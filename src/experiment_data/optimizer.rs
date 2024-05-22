use std::fs::{File, OpenOptions};
use std::io::Write;
use gkquad::prelude::Integrator;
use num::complex::Complex64;
use rand::{random, Rng, thread_rng};
use rand::rngs::ThreadRng;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator};
use rayon::iter::ParallelIterator;
use special::Gamma;
use crate::experiment_data::bootstrap_data::BootstrapData;


#[derive(Debug, Copy, Clone)]
pub struct Parameters<T> {
    pub a: T,
    pub b: T,
    pub n: T,
    pub d: T
}

impl Parameters<f64> {
    pub fn random_new() -> Self{
        let mut rng = rand::thread_rng();
        Parameters{
            a: rng.gen::<f64>() * 2.0 - 1.0,
            b: rng.gen::<f64>() * 6.0,
            n: rng.gen::<f64>() * 8.0 + 0.01,
            d: rng.gen::<f64>() * 4.0 - 1.0
        }
    }

    fn mutate_soft(x: f64, rng: &mut ThreadRng) -> f64 {
        x + (rng.gen::<f64>() * 0.1 - 0.05)
    }

    fn mutate_hard(x: f64, rng: &mut ThreadRng) -> f64 {
        x + (rng.gen::<f64>() * 1.0 - 0.5)
    }

    fn mutate_number(x: f64, rng: &mut ThreadRng) -> f64 {
        match rng.gen::<f64>() > 0.5 {
            true => Self::mutate_hard(x, rng),
            false => Self::mutate_soft(x, rng)
        }
    }

    pub fn mutate_random_param(&self, rng: &mut ThreadRng) -> Self {
        let num: f64 = rng.gen();
        match num {
            a if a > 0.75 => Parameters {
                a: Self::mutate_number(self.a, &mut thread_rng()),
                b: self.b,
                d: self.d,
                n: self.n
            },
            b if b > 0.5 => Parameters {
                a: self.a,
                b: Self::mutate_number(self.b, &mut thread_rng()),
                d: self.d,
                n: self.n
            },
            d if d > 0.25 => Parameters {
                a: self.a,
                b: self.b,
                d: Self::mutate_number(self.d, &mut thread_rng()),
                n: self.n
            },
            _ => Parameters {
                a: self.a,
                b: self.b,
                d: self.d,
                n: Self::mutate_number(self.n, &mut thread_rng())
            },
        }
    }

    pub fn mutate(&self) -> Self{
        let mut rng = rand::thread_rng();
        let n: u8 = rng.gen::<u8>() % 2 + 1;
        let mut params = self.clone();
        for i in 0..n {
            params = params.mutate_random_param(&mut rng);
        }
        params
    }

    pub fn merge(&self, other: &Self) -> Self {
        Parameters{
            a: (self.a + other.a) / 2.0,
            b: (self.b + other.b) / 2.0,
            n: (self.n + other.n) / 2.0,
            d: (self.d + other.d) / 2.0
        }
    }
}

pub struct Optimizer {
    q_data: Vec<(f64, BootstrapData, Complex64)>,
    population: Vec<Parameters<f64>>,
    real: bool,
    step: usize,
    last_best: (Parameters<f64>, f64)
}

impl Optimizer {

    pub fn new(q_data: Vec<(f64, BootstrapData)>, real: bool) -> Self{
        Optimizer{ q_data: q_data.into_iter().map(|(v, vec)|{
            let error = vec.boot_error_squared();
            (v, vec, error)
        }).collect(),
            real, population: (0..500).map(|_|Parameters::random_new()).collect(), step: 0, last_best: (Parameters{a:0.0,b:0.0,n:0.0,d:0.0}, 10e15)}
    }

    pub fn population_step(&mut self, sample_num: usize) -> f64{
        self.step += 1;
        let mut population_with_values: Vec<(Parameters<f64>, f64)> = self.population
            .par_iter()
            .map(|p|(*p, match self.real {
                true => self.calculate_qv(p.a, p.b, p.d, sample_num),
                false =>self.calculate_qv2s(p.a, p.b, p.n, p.d, sample_num)
            }))
            .collect();
        population_with_values.sort_by(|a,b|a.1.total_cmp(&b.1));
        let res = self.last_best.1 - population_with_values[0].1;
        self.last_best = (population_with_values[0].0, population_with_values[0].1);
        //println!("SAMPLE {} Best from this population {:?} with value: {} vs last value diff: {}", self.step, population_with_values[0].0, population_with_values[0].1, res);
        let mut to_append: Vec<Parameters<f64>> = Vec::new();
        to_append.push(population_with_values[0].0);
        let mut new_population= population_with_values.iter()
            .take(100)
            .map(|(k,_v)|k.mutate())
            .collect::<Vec<Parameters<f64>>>();
        let mut rng = thread_rng();
        for _ in 0..299 {
            let i: usize = rng.gen::<usize>() % 100;
            let j :usize = rng.gen::<usize>() % 100;
            to_append.push(population_with_values[i].0.merge(&population_with_values[j].0));
        }
        for _ in 0..100 {
            let i: usize = rng.gen::<usize>() % 100;
            let j :usize = rng.gen::<usize>() % 100;
            to_append.push(population_with_values[i].0.merge(&population_with_values[j].0).mutate());
        }
        to_append.append(&mut new_population);
        self.population = to_append;
        res
    }

    pub fn optimize(&mut self, max_steps: usize, delta: f64) -> Parameters<BootstrapData> {
        let size = self.q_data.get(0).unwrap().1.len();
        let mut i = 0;
        let mut vec_a = Vec::new();
        let mut vec_b = Vec::new();
        let mut vec_n = Vec::new();
        let mut vec_d = Vec::new();
        let mut vec_average = Vec::new();
        let mut file_av = OpenOptions::new().create(true).append(true).open("averages.dat").unwrap();
        while i < size {
            let res = self.optimize_one(max_steps, delta, i);
            vec_average.push(res.1);
            self.step = 0;
            vec_a.push(Complex64::new(res.0.a, 0.0));
            vec_b.push(Complex64::new(res.0.b, 0.0));
            vec_n.push(Complex64::new(res.0.n, 0.0));
            vec_d.push(Complex64::new(res.0.d, 0.0));
            println!("Optimized index {}, with result {:?}", i, res);
            i+=1;
        }
        let res = self.optimize_one(max_steps, delta, i);
        file_av.write(format!("Average error: {}", vec_average.iter().sum::<f64>() / vec_average.len() as f64).as_bytes()).expect("TODO: panic message");
        println!("\nAverage error: {}\n", vec_average.iter().sum::<f64>() / vec_average.len() as f64);
        Parameters{
            a: BootstrapData::new(vec_a, Complex64::new(res.0.a, 0.0)),
            b: BootstrapData::new(vec_b, Complex64::new(res.0.b, 0.0)),
            n: BootstrapData::new(vec_n, Complex64::new(res.0.n, 0.0)),
            d: BootstrapData::new(vec_d, Complex64::new(res.0.d, 0.0)),
        }
    }

    pub fn optimize_one(&mut self, max_steps: usize, delta: f64, index:usize) -> (Parameters<f64>, f64) {
        let mut last_jump = 0;
        self.population = (0..500).map(|_|Parameters::random_new()).collect();
        loop {
            let current_delta = self.population_step(index);
            if current_delta > delta{
                last_jump = 0
            } else {
                last_jump += 1;
            }
            let step = self.step;
            if step >= max_steps || last_jump >= 20 {
                break;
            }
        }
        self.last_best
    }
    pub fn integral_a_b_n_cos(a: f64, b:f64, d:f64, n: f64, ni: f64) -> f64{
        unsafe {
            let v = Integrator::new(|x: f64|{
                (x*ni).cos() * n * x.powf(a) * (1.0-x).powf(b) * (1.0 + d * x.sqrt())
            }).max_iters(3).run(0.0..1.0).estimate_unchecked();
            v
        }
    }

    pub fn integral_a_b_n_sin(a: f64, b:f64, n: f64, d: f64, ni: f64) -> f64{
        unsafe {
            let v = Integrator::new(|x: f64|{
                (x*ni).sin() * n * x.powf(a) * (1.0-x).powf(b) * (1.0 + d * x.sqrt())
            }).max_iters(3).run(0.0..1.0).estimate_unchecked();
            v
        }
    }
    pub fn calculate_chi_elem_a_b_n(ni: f64, data: &BootstrapData, a: f64, b:f64, n: f64, d:f64, real: bool, sample_num: usize, error: Complex64) -> f64 {
        match real {
            true => {
                let int_res = Self::integral_a_b_n_cos(a,b,n,d,ni);
                let res = (data.get_sample(sample_num).re - int_res).powi(2)/error.re;
                res
            }
            false => {
                let int_res = Self::integral_a_b_n_sin(a,b,n,d,ni);
                let res = (data.get_sample(sample_num).im - int_res).powi(2)/error.im;
                res
            }
        }
    }

    pub fn calculate_chi_qv(data: &Vec<(f64, BootstrapData, Complex64)>, a: f64, b:f64, d:f64, sample_num: usize) -> f64 {
        let c_part =(Gamma::gamma(a+1.5)*Gamma::gamma(b+1.0)) * d / Gamma::gamma(a + b + 2.5);
        let normal_part = (Gamma::gamma(a+1.0)*Gamma::gamma(b+1.0))/ Gamma::gamma(a + b + 2.0);
        let n = 1.0/(c_part + normal_part);
        let res = data
            .iter()
            .map(|(ni,data, error)|Self::calculate_chi_elem_a_b_n(*ni, data, a, b, n,d, true,sample_num, *error))
            .sum();
        res
    }

    pub fn calculate_chi_qv2s(data: &Vec<(f64, BootstrapData, Complex64)>, a: f64, b:f64, d: f64, n:f64, sample_num: usize) -> f64 {
        let res = data
            .iter()
            .map(|(ni,data, error)|Self::calculate_chi_elem_a_b_n(*ni, data, a, b, n,d, false,sample_num, *error))
            .sum();
        res
    }

    pub fn calculate_qv(&self, a: f64, b:f64, d: f64, sample_num: usize) -> f64 {
        Self::calculate_chi_qv(&self.q_data, a, b, d, sample_num)
    }

    pub fn calculate_qv2s(&self, a: f64, b:f64, n: f64, d: f64, sample_num: usize) -> f64 {
        Self::calculate_chi_qv2s(&self.q_data, a, b, n, d, sample_num)
    }
}