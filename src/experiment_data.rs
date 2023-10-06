mod bootstrap_data;

use std::collections::HashMap;
use std::f64::consts::PI;
use std::{fs, iter};
use std::fs::{DirEntry, File};
use std::path::Path;
use num::complex::Complex64;
use rayon::iter::IntoParallelRefIterator;
use crate::experiment_data::bootstrap_data::BootstrapData;

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
    data: HashMap<DataInfo, BootstrapData>,
    reduced_m: HashMap<DataInfo, BootstrapData>
    //fits: Option<HashMap<u8, Vec<BootstrapData>>>
}

impl ExperimentData {

    const LATTICE_LENGTH: i32 = 48;

    pub fn load_and_calculate(data_path: &Path) -> ExperimentData {
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
        ExperimentData{ max_mom, max_z, data, reduced_m }
    }

    pub fn convert_to_ioffe_time(mom: u8, z:u8) -> f64 {
        mom as f64 * z as f64 * 2.0 * PI / Self::LATTICE_LENGTH as f64
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

    fn calculate_reduced_m_value(data: &HashMap<DataInfo, BootstrapData>, mom: u8, z:u8) -> BootstrapData{
        let bare = Self::get_from(data, mom, z);
        let bare_0_z = Self::get_from(data, mom, 0);
        let bare_0_v = Self::get_from(data, 0, z);
        let bare_0_0 = Self::get_from(data, 0, 0);
        ((bare / bare_0_z) / (bare_0_v / bare_0_0))
    }

    fn calculate_reduced_m(data: &HashMap<DataInfo, BootstrapData>) -> HashMap<DataInfo, BootstrapData>{
        data.iter()
            .map(|(k,v)|(*k, Self::calculate_reduced_m_value(data, k.mom, k.z)))
            .collect()
    }

    fn get_from(data: &HashMap<DataInfo, BootstrapData>, mom: u8, z:u8) -> &BootstrapData {
        data.get(&DataInfo::new(mom, z)).unwrap()
    }

    /*pub fn fit(&self, z: u8, order: usize) -> Vec<BootstrapData> {
        let order = order * 2 - 1;
        let x_values: Vec<f64> = (0..self.max_mom +1).map(|e|e as f64).collect();
        let full_reduced = self.get_full_reduced_m();
        let y_values: Vec<&BootstrapData> = (0..self.max_mom +1).map(|i|Self::get_from(&full_reduced, i,z)).collect();
        BootstrapData::perform_operation_multiple_to_multiple(y_values, |vec|{
            let (y_re, y_im): (Vec<f64>, Vec<f64>) = vec.iter().map(|e|(e.re, e.im)).unzip();
            let re_res = polyfit(&x_values, &y_re, order).expect("Can not fit real part");
            let im_res = polyfit(&x_values, &y_im, order).expect("Can not fit real part");
            zip(re_res, im_res).map(|(a,b)|Complex64::new(a,b)).collect()
        })
    }*/


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