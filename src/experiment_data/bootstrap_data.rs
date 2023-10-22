use std::ops::{Add, Div, Mul, Sub};
use num::complex::{Complex64};
use rayon::prelude::*;
#[derive(Debug, Clone)]
pub struct BootstrapData {
    sample_0: Complex64,
    bootstrap_samples: Vec<Complex64>
}

impl BootstrapData {

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

    pub fn value(&self) -> Complex64 {
        self.sample_0
    }

    pub fn value_error(&self) -> (Complex64, Complex64) {
        (self.sample_0, self.boot_error())
    }

    pub fn value_average_error(&self) -> (Complex64, Complex64, Complex64) {
        (self.sample_0, self.boot_average(), self.boot_error())
    }

    pub fn perform_operation<F>(&self, operation: F) -> BootstrapData
        where F: Fn(Complex64) -> Complex64 {
        BootstrapData {
            sample_0: operation(self.sample_0),
            bootstrap_samples: self.bootstrap_samples.iter().map(|e|operation(*e)).collect()
        }
    }

    pub fn perform_operation_with_other<F>(&self, other: &BootstrapData, operation: F) -> BootstrapData
        where F: Fn(Complex64, Complex64) -> Complex64 {
        BootstrapData {
            sample_0: operation(self.sample_0, other.sample_0),
            bootstrap_samples: self.bootstrap_samples.iter().zip(other.bootstrap_samples.iter()).map(|(a,b)|operation(*a,*b)).collect()
        }
    }

    pub fn perform_operation_multiple<F>(data: Vec<&BootstrapData>, operation: F) -> BootstrapData
        where F: Fn(Vec<Complex64>) -> Complex64 {
        BootstrapData {
            sample_0: operation(data.iter().map(|e|e.sample_0).collect()),
            bootstrap_samples: (0..data[0].bootstrap_samples.len()).map(|i|{
                operation(data.iter().map(|e|e.bootstrap_samples[i]).collect())
            }).collect()
        }
    }

    pub fn perform_operation_multiple_to_multiple<F>(data: Vec<&BootstrapData>, operation: F) -> Vec<BootstrapData>
        where F: Fn(Vec<Complex64>) -> Vec<Complex64> {
        let samples_0 = operation(data.iter().map(|e|e.sample_0).collect());
        let bootstrap_samples: Vec<Vec<Complex64>> = (0..data[0].bootstrap_samples.len()).map(|i|{
            operation(data.iter().map(|e|e.bootstrap_samples[i]).collect())
        }).collect();
        (0..samples_0.len())
            .map(|i|BootstrapData{
                sample_0: samples_0[i],
                bootstrap_samples: bootstrap_samples[i].clone()
            })
            .collect()
    }

}

impl From<Vec<Complex64>> for BootstrapData{
    fn from(value: Vec<Complex64>) -> Self {
        let mut value = value.into_iter();
        let sample_0 = value.next().expect("Empty vector while creating bootstrap data");
        let bootstrap_samples = value.collect();
        BootstrapData { sample_0, bootstrap_samples }
    }
}

impl Div for &BootstrapData {
    type Output = BootstrapData;

    fn div(self, rhs: Self) -> Self::Output {
        self.perform_operation_with_other(rhs, |a,b|a/b)
    }
}

impl Div for BootstrapData {
    type Output = BootstrapData;

    fn div(self, rhs: Self) -> Self::Output {
        self.perform_operation_with_other(&rhs, |a,b|a/b)
    }
}

impl Sub for &BootstrapData {
    type Output = BootstrapData;
    fn sub(self, rhs: Self) -> Self::Output { self.perform_operation_with_other(rhs, |a,b|a-b) }
}

impl Sub for BootstrapData {
    type Output = BootstrapData;
    fn sub(self, rhs: Self) -> Self::Output { self.perform_operation_with_other(&rhs, |a,b|a-b) }
}

impl Mul<f64> for BootstrapData {
    type Output = BootstrapData;
    fn mul(self, rhs: f64) -> Self::Output { self.perform_operation(|a|a*rhs) }
}

impl Add for &BootstrapData {
    type Output = BootstrapData;

    fn add(self, rhs: Self) -> Self::Output {
        self.perform_operation_with_other(rhs, |a,b|a+b)
    }
}

impl Add for BootstrapData {
    type Output = BootstrapData;

    fn add(self, rhs: Self) -> Self::Output {
        self.perform_operation_with_other(&rhs, |a,b|a+b)
    }
}