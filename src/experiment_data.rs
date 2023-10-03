mod bootstrap_data;
mod bootstrap_data_creator;

use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs;
use std::fs::{DirEntry, File};
use std::iter::Map;
use std::path::Path;
use num::complex::Complex64;
use crate::experiment_data::bootstrap_data::BootstrapData;

#[derive(PartialEq, Eq, Hash, Clone)]
enum ZDirection{
    Plus, Minus, Average
}
#[derive(PartialEq, Eq, Hash, Clone)]
struct DataInfo{
    z_direction: ZDirection,
    mom: u8,
    z: u8
}

pub struct ExperimentData{
    max_mom: u8,
    max_z: u8,
    data: HashMap<DataInfo, BootstrapData>
}

impl ExperimentData {

    pub fn load(data_path: &Path) -> ExperimentData {
        let max_mom = 0u8;
        let max_z = 0u8;
        data_path
            .read_dir()
            .expect(format!("Load experiment data failed '{}' is not a dir or can not be read", data_path.display()).as_str())
            .filter_map(Result::ok)
            .map(ExperimentData::read_file)

    }

    fn read_file(file: DirEntry) -> Map<DataInfo, BootstrapData> {
        let file_input = fs::read_to_string(file).expect(format!("Can not read file '{}' is not a dir or can not be read", file.path().display()).as_str());
        let (z, mom) = Self::read_z_and_mom_from_filename(file.file_name().as_os_str().to_str().expect(&format!("File name {} incorrect", file.path().display())));
        let mut complex_stream:Vec<(Complex64, Complex64, Complex64)> = file_input.lines()
            .into_iter()
            .map(ExperimentData::input_line_to_complex_numbers)
            .collect();

        let sample_0 = complex_stream.next().unwrap();
        let bootstrap_samples = complex_stream.collect();
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