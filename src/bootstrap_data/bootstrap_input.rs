use num::complex::{Complex64};
use std::{fs};
use std::error::Error;
use std::string::ParseError;
use rayon::prelude::*;
#[derive(Debug)]
pub struct BootstrapInput {
    sample_0: Complex64,
    bootstrap_samples: Vec<Complex64>
}

impl BootstrapInput {
    pub fn read_from_file(file_path: &str, plus_sign: bool) -> Result<BootstrapInput, Box<dyn Error>>{
        let file_input = fs::read_to_string(file_path)?;
        let mut complex_stream = file_input.lines()
            .into_iter()
            .map(|e|Self::input_line_to_complex_num(e, plus_sign))
            .map(|e|e.unwrap());
        let sample_0 = complex_stream.next().unwrap();
        let bootstrap_samples = complex_stream.collect();
        Ok(BootstrapInput{ sample_0, bootstrap_samples })
    }

    pub fn boot_average(&self) -> Complex64 {
        self.bootstrap_samples.par_iter().sum::<Complex64>() / (self.bootstrap_samples.len() as f64)
    }

    pub fn get_sample_0(&self) -> Complex64 {
        self.sample_0
    }

    pub fn boot_error(&self) -> Complex64 {
        let boot_average = self.boot_average();
        let sample_size: f64 = self.bootstrap_samples.len() as f64;
        let res_sum: (f64,f64) = self.bootstrap_samples
            .par_iter()
            .map(|e|(
                (e.re - boot_average.re).powi(2)/sample_size,
                (e.im - boot_average.im).powi(2)/sample_size
            ))
            .reduce(||{(0.0,0.0)}, |(a,b), (c,d)|{
                (a+c, b+d)
            });
        Complex64::new(res_sum.0.sqrt(), res_sum.1.sqrt())
    }
    fn input_line_to_complex_num(line: &str, plus_sign : bool) -> Result<Complex64, ParseError>{
        let mut line_split = line.split_whitespace();
        let num = match plus_sign {
            true => 4,
            false => 2
        };
        let mut line_split = line_split.skip(num);
        let re = line_split.next().map(|s|s.parse::<f64>().unwrap()).unwrap();
        let im = line_split.next().map(|s|s.parse::<f64>().unwrap()).unwrap();
        Ok(Complex64::new(re, im))
    }

}