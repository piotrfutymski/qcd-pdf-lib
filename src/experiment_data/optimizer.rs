use std::cmp::max;
use std::collections::HashMap;
use gkquad::prelude::Integrator;
use rand::{Rng, thread_rng};
use rand::rngs::ThreadRng;
use special::Gamma;
use crate::experiment_data::bootstrap_data::BootstrapData;


#[derive(Debug, Copy, Clone)]
pub struct Parameters {
    a: f64,
    b: f64,
    n: f64
}

impl Parameters {
    pub fn random_new() -> Self{
        let mut rng = rand::thread_rng();
        Parameters{
            a: rng.gen::<f64>() * 2.9 + 0.1,
            b: rng.gen::<f64>() * 2.9 + 0.1,
            n: rng.gen::<f64>() * 2.9 + 0.1
        }
    }

    fn mutate_soft(x: f64, rng: &mut ThreadRng) -> f64 {
        f64::max(rng.gen::<f64>() * 0.01 - 0.005 + x, 0.001)
    }

    fn mutate_hard(x: f64, rng: &mut ThreadRng) -> f64 {
        f64::max(rng.gen::<f64>() * 1.0 - 0.5 + x, 0.001)
    }

    fn mutate_number(x: f64, rng: &mut ThreadRng) -> f64 {
        match rng.gen::<f64>() > 0.8 {
            true => Self::mutate_hard(x, rng),
            false => Self::mutate_soft(x, rng)
        }
    }

    pub fn mutate(&self) -> Self{
        let mut rng = rand::thread_rng();
        Parameters{
            a: Self::mutate_number(self.a, &mut rng),
            b: Self::mutate_number(self.b, &mut rng),
            n: Self::mutate_number(self.n, &mut rng),
        }
    }

    pub fn merge(&self, other: &Self) -> Self {
        Parameters{
            a: (self.a + other.a) / 2.0,
            b: (self.b + other.b) / 2.0,
            n: (self.n + other.n) / 2.0
        }
    }
}

pub struct Optimizer {
    q_data: Vec<(f64, BootstrapData)>,
    population: Vec<Parameters>,
    real: bool,
    step: usize,
    last_best: (Parameters, f64)
}

impl Optimizer {

    pub fn new(q_data: Vec<(f64, BootstrapData)>, real: bool) -> Self{
        Optimizer{q_data, real, population: (0..200).map(|_|Parameters::random_new()).collect(), step: 0, last_best: (Parameters{a:0.0,b:0.0,n:0.0}, 10e15)}
    }

    pub fn population_step(&mut self) -> f64{
        self.step += 1;
        println!("Calculating best population for step {}", self.step);
        let eval_func:Box<dyn Fn(&Parameters) -> f64> = match self.real {
            true => Box::new(|p:&Parameters| self.calculate_qv(p.a, p.b)),
            false => Box::new(|p:&Parameters| self.calculate_qv2s(p.a, p.b, p.n))
        };
        let mut population_with_values: Vec<(Parameters, f64)> = self.population
            .iter()
            .map(|e|(*e, eval_func(e)))
            .collect();
        population_with_values.sort_by(|a,b|a.1.total_cmp(&b.1));
        drop(eval_func);
        let res = self.last_best.1 - population_with_values[0].1;
        self.last_best = (population_with_values[0].0, population_with_values[0].1);
        println!("Best from this population {:?} with value: {} vs last value diff: {}", population_with_values[0].0, population_with_values[0].1, res);
        println!("Mutating best 20 parameter sets");
        let mut to_append: Vec<Parameters> = Vec::new();
        to_append.push(population_with_values[0].0);
        let mut new_population= population_with_values.iter()
            .take(40)
            .map(|(k,v)|k.mutate())
            .collect::<Vec<Parameters>>();
        println!("Merging randomly for new population");
        let mut rng = thread_rng();
        for _ in 0..159 {
            let i: usize = rng.gen::<usize>() % 40;
            let j :usize = rng.gen::<usize>() % 40;
            to_append.push(new_population[i].merge(&new_population[j]));
        }
        to_append.append(&mut new_population);
        self.population = to_append;
        res
    }

    pub fn optimize(&mut self, max_steps: usize, delta: f64) -> Parameters {
        let mut deltas = [10e15, 10e15, 10e15];
        loop {
            let current_delta = self.population_step();
            if current_delta < deltas[2] {
                deltas = [deltas[0], deltas[1], current_delta];
                deltas.sort_by(|a,b|a.total_cmp(b));
            }
            let step = self.step;
            if step >= max_steps {
                break;
            }
        }
        self.last_best.0
    }
    pub fn integral_a_b_n_cos(a: f64, b:f64, n: f64, ni: f64) -> f64{
        unsafe {
            let v = Integrator::new(|x: f64|{
                (x*ni).cos() * n * x.powf(a) * (1.0-x).powf(b)
            }).run(0.0..1.0).estimate_unchecked();
            v
        }
    }

    pub fn integral_a_b_n_sin(a: f64, b:f64, n: f64, ni: f64) -> f64{
        unsafe {
            let v = Integrator::new(|x: f64|{
                (x*ni).sin() * n * x.powf(a) * (1.0-x).powf(b)
            }).run(0.0..1.0).estimate_unchecked();
            v
        }
    }
    pub fn calculate_chi_elem_a_b_n(ni: f64, data: &BootstrapData, a: f64, b:f64, n: f64, real: bool) -> f64 {
        match real {
            true => {
                let int_res = Self::integral_a_b_n_cos(a,b,n,ni);
                let res = (data.boot_average().re - int_res).powi(2)/data.boot_error_squared().re;
                res
            }
            false => {
                let int_res = Self::integral_a_b_n_sin(a,b,n,ni);
                let res = (data.boot_average().im - int_res).powi(2)/data.boot_error_squared().im;
                res
            }
        }
    }

    pub fn calculate_chi_qv(data: &Vec<(f64, BootstrapData)>, a: f64, b:f64) -> f64 {
        let n = Gamma::gamma(a + b + 2.0) / (Gamma::gamma(a+1.0)*Gamma::gamma(b+1.0));
        let res = data
            .iter()
            .map(|(ni,d)|Self::calculate_chi_elem_a_b_n(*ni, d, a, b, n, true))
            .sum();
        res
    }

    pub fn calculate_chi_qv2s(data: &Vec<(f64, BootstrapData)>, a: f64, b:f64, n:f64) -> f64 {
        let res = data
            .iter()
            .map(|(ni,d)|Self::calculate_chi_elem_a_b_n(*ni, d, a, b, n, false))
            .sum();
        res
    }

    pub fn calculate_qv(&self, a: f64, b:f64) -> f64 {
        Self::calculate_chi_qv(&self.q_data, a, b)
    }

    pub fn calculate_qv2s(&self, a: f64, b:f64, n: f64) -> f64 {
        Self::calculate_chi_qv2s(&self.q_data, a, b, n)
    }
}