mod bootstrap_data;
mod optimizer;

use std::collections::HashMap;
use std::f64::consts::PI;
use std::{fs, iter};
use std::fs::{DirEntry, File};
use std::hint::black_box;
use std::io::Write;
use std::path::Path;
use gkquad::single::Integrator;
use num::complex::{Complex64, ComplexFloat};
use rayon::iter::IntoParallelRefIterator;
use crate::experiment_data::bootstrap_data::BootstrapData;
use crate::experiment_data::optimizer::{Optimizer, Parameters};

#[derive(PartialEq, Eq, Hash, Clone, Default, Copy)]
pub enum ZDirection{
    Plus,
    Minus,
    #[default]
    Average
}
#[derive(PartialEq, Eq, Hash, Clone, Copy)]
pub struct DataInfo{
    z_direction: ZDirection,
    mom: u8,
    z: u8
}

impl DataInfo {
    pub fn new(mom: u8, z:u8) -> DataInfo{
        DataInfo{ z_direction: ZDirection::default(), mom, z}
    }
}

pub struct ExperimentData{
    max_mom: u8,
    max_z: u8,
    max_num_to_average: u8,
    data: HashMap<DataInfo, BootstrapData>,
    reduced_m: HashMap<DataInfo, BootstrapData>,
    m_prime: HashMap<DataInfo, BootstrapData>,
    q: HashMap<DataInfo, BootstrapData>,
    q_averaged: HashMap<u8, BootstrapData>,
    best_params_qv: Parameters,
    best_params_qv2s: Parameters
}

impl ExperimentData {

    const LATTICE_LENGTH: i32 = 48;
    const ALFA_S_PI: f64 = 0.1287;
    const C_F: f64 = 4.0 / 3.0;

    const MI_A: f64 = 0.950711;

    const GAMMA_E: f64 = 0.5772156649;

    pub fn load_and_calculate(data_path: &Path, max_num_to_average: u8) -> ExperimentData {
        let data: HashMap<DataInfo, BootstrapData> = data_path
            .read_dir()
            .expect(format!("Load experiment data failed '{}' is not a dir or can not be read", data_path.display()).as_str())
            .filter_map(Result::ok)
            .map(ExperimentData::read_file)
            .flatten()
            .collect();
        let max_mom = data.keys().map(|e|e.mom).max().unwrap_or_default();
        let max_z = data.keys().map(|e|e.z).max().unwrap_or_default();
        let reduced_m = Self::calculate_reduced_m(&data);
        let m_prime = Self::calculate_m_prime(max_mom, &reduced_m);
        let q = Self::calculate_q(max_mom, &m_prime, &reduced_m);
        let q_averaged = Self::calculate_q_averaged(&q, max_num_to_average);
        let best_params_qv = Self::calculate_best_params(&q_averaged, true);
        let best_params_qv2s = Self::calculate_best_params(&q_averaged, false);
        ExperimentData{ max_mom, max_z, data, reduced_m, m_prime, q, max_num_to_average, q_averaged, best_params_qv, best_params_qv2s}
    }

    pub fn convert_to_ioffe_time(mom: f64, z:u8) -> f64 {
        mom * z as f64 * 2.0 * PI / Self::LATTICE_LENGTH as f64
    }

    pub fn convert_to_ioffe_time_one_arg(ni: u8) -> f64 {
        ni as f64 * 2.0 * PI / Self::LATTICE_LENGTH as f64
    }

    pub fn convert_from_ioffe_time(it: f64, z:u8) -> f64 {
        it * Self::LATTICE_LENGTH as f64 / (2.0 * PI * z as f64)
    }

    pub fn bare_m(&self) -> &HashMap<DataInfo, BootstrapData>{
        &self.data
    }

    pub fn bare_m_value(&self, mom: u8, z:u8) -> &BootstrapData {
        Self::get_from(self.bare_m(), mom, z)
    }

    pub fn reduced_m(&self) -> &HashMap<DataInfo, BootstrapData> {
        &self.reduced_m
    }

    pub fn reduced_m_value(&self, mom: u8, z:u8) -> &BootstrapData {
        Self::get_from(self.reduced_m(), mom, z)
    }

    pub fn m_prime(&self) -> &HashMap<DataInfo, BootstrapData> {
        &self.m_prime
    }

    pub fn m_prime_value(&self, mom: u8, z:u8) -> &BootstrapData {
        Self::get_from(self.m_prime(), mom, z)
    }

    pub fn q(&self) -> &HashMap<DataInfo, BootstrapData> {
        &self.q
    }

    pub fn q_value(&self, mom: u8, z:u8) -> &BootstrapData {
        Self::get_from(self.q(), mom, z)
    }


    pub fn q_averaged(&self) -> &HashMap<u8, BootstrapData> {
        &self.q_averaged
    }
    pub fn q_averaged_value(&self, ni: u8) -> &BootstrapData {
        self.q_averaged.get(&ni).unwrap()
    }

    pub fn m_prime_value_it(&self, it: f64, z:u8) -> &BootstrapData {
        let mom = Self::convert_from_ioffe_time(it, z).round() as u8;
        Self::get_from(self.m_prime(), mom, z)
    }

    pub fn generate_output_plot_file(data: &HashMap<DataInfo, BootstrapData>, mom:u8, sample_count:u8, filename:&str, re: bool) {
        let mut file = File::create(filename).unwrap();
        for z in 0..sample_count + 1 {
            let x = Self::convert_to_ioffe_time(mom as f64, z);
            let bootstrap_data = Self::get_from(data, mom, z);
            let (y, err) = match re {
                true => (bootstrap_data.boot_average().re, bootstrap_data.boot_error().re),
                false => (bootstrap_data.boot_average().im, bootstrap_data.boot_error().im)
            };
            file.write(format!("{} {} {}\n", x, y, err).as_bytes()).unwrap();
        }
    }

    pub fn generate_output_plot_file_one_arg(data: &HashMap<u8, BootstrapData>, sample_count:u8, filename:&str, re: bool) {
        let mut file = File::create(filename).unwrap();
        for ni in 0..sample_count + 1 {
            let x = Self::convert_to_ioffe_time_one_arg(ni);
            if let Some(bootstrap_data) = data.get(&ni) {
                let (y, err) = match re {
                    true => (bootstrap_data.boot_average().re, bootstrap_data.boot_error().re),
                    false => (bootstrap_data.boot_average().im, bootstrap_data.boot_error().im)
                };
                file.write(format!("{} {} {}\n", x, y, err).as_bytes()).unwrap();
            }
        }
    }

    fn calculate_reduced_m_value(data: &HashMap<DataInfo, BootstrapData>, mom: u8, z:u8) -> BootstrapData{
        let bare = Self::get_from(data, mom, z);
        let bare_0_z = Self::get_from(data, mom, 0);
        let bare_0_v = Self::get_from(data, 0, z);
        let bare_0_0 = Self::get_from(data, 0, 0);
        (bare / bare_0_z) / (bare_0_v / bare_0_0)
    }

    fn calculate_reduced_m(data: &HashMap<DataInfo, BootstrapData>) -> HashMap<DataInfo, BootstrapData>{
        data.iter()
            .filter(|(k,v)|k.z_direction == ZDirection::Average)
            .map(|(k,v)|(*k, Self::calculate_reduced_m_value(data, k.mom, k.z)))
            .collect()
    }

    fn calculate_m_prime_value(max_mom: u8, data: &HashMap<DataInfo, BootstrapData>, mom: u8, z:u8) -> BootstrapData {
        println!("Calculating m prime value for mom={}, z={}", mom, z);
        let vec: Vec<&BootstrapData> = (0..max_mom+1)
            .map(|e|Self::get_from(data, e, z))
            .collect();
        BootstrapData::perform_operation_multiple(&vec, |c_data|{ Self::integral_m_prim(&c_data, z, mom)})
    }

    fn calculate_q_value(max_mom: u8, m_prime: &HashMap<DataInfo, BootstrapData>, m: &HashMap<DataInfo, BootstrapData>, mom: u8, z:u8) -> BootstrapData {
        println!("Calculating q value for mom={}, z={}", mom, z);
        let mut vec: Vec<&BootstrapData> = (0..max_mom+1)
            .map(|e|Self::get_from(m, e, z))
            .collect();
        vec.push(Self::get_from(m_prime, mom, z));
        BootstrapData::perform_operation_multiple(&vec, |c_data|{ Self::integral_q(&c_data, z, mom)})
    }

    fn calculate_q_averaged_value(data: &Vec<&BootstrapData>) -> BootstrapData {
        BootstrapData::perform_operation_multiple(data, |c_data|{
            match c_data.len() {
                1 => c_data[0],
                _ => c_data.iter().sum::<Complex64>() / (c_data.len() as f64)
            }
        })
    }
    fn calculate_m_prime(max_mom: u8, data: &HashMap<DataInfo, BootstrapData>) -> HashMap<DataInfo, BootstrapData>{
        data.iter()
            .filter(|(k,v)|k.z_direction == ZDirection::Average)
            .map(|(k,v)|(*k, Self::calculate_m_prime_value(max_mom, data, k.mom, k.z)))
            .collect()
    }

    fn calculate_q(max_mom: u8, m_prime: &HashMap<DataInfo, BootstrapData>, m: &HashMap<DataInfo, BootstrapData>) -> HashMap<DataInfo, BootstrapData>{
        m_prime.iter()
            .filter(|(k,v)|k.z_direction == ZDirection::Average)
            .map(|(k,v)|(*k, Self::calculate_q_value(max_mom, m_prime, m, k.mom, k.z)))
            .collect()
    }

    fn calculate_q_averaged(q: &HashMap<DataInfo, BootstrapData>, max_z: u8) -> HashMap<u8, BootstrapData> {
        let mut grouped_map = HashMap::new();

        q.iter()
            .filter(|(k,v)|k.z_direction == ZDirection::Average && k.z <= max_z && k.mom != 0)
            .for_each(|(k,v)|{
                let num = k.z *k.mom;
                grouped_map.entry(num).or_insert(Vec::new()).push(v)
            });

        grouped_map
            .iter()
            .map(|(k,v)|(*k, Self::calculate_q_averaged_value(v)))
            .collect()
    }

    fn calculate_best_params(q_averaged: &HashMap<u8, BootstrapData>, real: bool) -> Parameters {
        let mut optimizer = Optimizer::new(
            q_averaged.iter()
                .filter(|(k,v)|**k != 0)
                .map(|(k,v)|(Self::convert_to_ioffe_time_one_arg(*k), v.clone()))
                .collect(),
            real
        );
        optimizer.optimize(100,0.01)
    }

    fn get_from(data: &HashMap<DataInfo, BootstrapData>, mom: u8, z:u8) -> &BootstrapData {
        data.get(&DataInfo::new(mom, z)).unwrap()
    }

    fn linear_interpolation(data: &Vec<Complex64>, v: f64) -> Complex64 {
        let min_v = v.floor() as u8;
        let max_v = min_v + 1;
        let frac = v.fract();
        if max_v >= data.len() as u8 {
            return *data.last().unwrap();
        }
        let up_v: Complex64 = data[max_v as usize];
        let down_v: Complex64 = data[min_v as usize];
        down_v + ((up_v - down_v) * frac)
    }

    fn integral_m_prim(data: &Vec<Complex64>, z: u8, mom: u8) -> Complex64 {
        let m: Complex64 = data[mom as usize];
        let int_mul = - Self::C_F * Self::ALFA_S_PI / 2.0;
        let g_e = (2.0*Self::GAMMA_E + 1.0).exp();
        if z == 0 {
            return m;
        }
        let v = 0.25 * (z as f64).powi(2) * Self::MI_A.powi(2)* g_e;
        let v = v.ln();
        unsafe {
            let re = Integrator::new(|x: f64|{
                let b = (1.0+(x.powi(2)))/(1.0-x);
                let m_u = Self::linear_interpolation(data, x*mom as f64).re;
                let res = -v*b * (m_u - m.re);
                res
            }).run(0.0..1.0).estimate_unchecked();
            let im = Integrator::new(|x: f64|{
                let b = (1.0+(x.powi(2)))/(1.0-x);
                let m_u = Self::linear_interpolation(data, x*mom as f64).im;
                let res = -v*b * (m_u - m.im);
                res
            }).run(0.0..1.0).estimate_unchecked();
            let v = m + int_mul * Complex64::new(re, im);
            v
        }
    }

    fn integral_q(data: &Vec<Complex64>, z: u8, mom: u8) -> Complex64 {
        let mut data_cloned = data.clone();
        let m_prime = data_cloned.pop().unwrap();
        let m: Complex64 = data_cloned[mom as usize];
        let int_mul = - Self::C_F * Self::ALFA_S_PI / 2.0;
        unsafe {
            let re = Integrator::new(|x: f64|{
                let l = 4.0*((1.0-x).ln()/(1.0-x)) - 2.0*(1.0-x);
                let m_u = Self::linear_interpolation(data, x*mom as f64).re;
                let res = -l * (m_u - m.re);
                res
            }).run(0.0..1.0).estimate_unchecked();
            let im = Integrator::new(|x: f64|{
                let l = 4.0*((1.0-x).ln()/(1.0-x)) - 2.0*(1.0-x);
                let m_u = Self::linear_interpolation(data, x*mom as f64).im;
                let res = -l * (m_u - m.im);
                res
            }).run(0.0..1.0).estimate_unchecked();
            let v = m_prime + int_mul * Complex64::new(re, im);
            v
        }
    }

    fn read_file(file: DirEntry) -> HashMap<DataInfo, BootstrapData> {
        let file_input = fs::read_to_string(file.path()).expect(format!("Can not read file '{}' is not a dir or can not be read", file.path().display()).as_str());
        let (z, mom) = Self::read_z_and_mom_from_filename(file.file_name().as_os_str().to_str().expect(&format!("File name {} incorrect", file.path().display())));
        let all_numbers: Vec<(Complex64, Complex64, Complex64)> = file_input.lines()
            .into_iter()
            .map(ExperimentData::input_line_to_complex_numbers)
            .collect();
        let ( mut vec_avg, mut vec_min, mut vec_plus)  = (vec![], vec![], vec![]);
        for complex_tuple in all_numbers.into_iter(){
            vec_avg.push(complex_tuple.0);
            vec_min.push(complex_tuple.1);
            vec_plus.push(complex_tuple.2);
        }
        let mut res = HashMap::new();
        res.insert(DataInfo{ z_direction: ZDirection::Average, mom, z, }, vec_avg.into());
        res.insert(DataInfo{ z_direction: ZDirection::Minus, mom, z, }, vec_min.into());
        res.insert(DataInfo{ z_direction: ZDirection::Plus, mom, z, }, vec_plus.into());
        res
    }


    fn read_z_and_mom_from_filename(filename: &str) -> (u8, u8) {
        let parts: Vec<&str> = filename.split("_").collect();
        let z_part = parts[1];
        let mom_part = parts[3];
        (
            z_part.strip_prefix("z").unwrap().parse().unwrap(),
            mom_part.strip_prefix("mom").unwrap().parse().unwrap()
        )
    }

    fn input_line_to_complex_numbers(line: &str) -> (Complex64, Complex64, Complex64){
        let numbers: Vec<f64> = line
            .split_whitespace()
            .map(|s|s.parse::<f64>().unwrap())
            .collect();
        (
            Complex64::new(numbers[0], numbers[1]),
            Complex64::new(numbers[2], numbers[3]),
            Complex64::new(numbers[4], numbers[5])
        )
    }

}