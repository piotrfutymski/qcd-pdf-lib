use num::complex::{Complex64, ComplexFloat};

pub trait FunctionEstimator<'a> {
    fn new(data: &'a Vec<Complex64>) -> Self where Self: Sized;
    fn estimate(&self, x: f64) -> Complex64;

}

pub trait ToFunctionEstimator {
    fn build_estimator<'a>(&self, data: &'a Vec<Complex64>) -> Box<dyn FunctionEstimator<'a> + 'a>;
}

pub struct LinearInterpolationEstimator<'a>{
    data: &'a Vec<Complex64>
}

impl <'a> FunctionEstimator<'a> for LinearInterpolationEstimator<'a> {
    fn new(data: &'a Vec<Complex64>) -> Self {
        return Self{data}
    }

    fn estimate(&self, x: f64) -> Complex64 {
        let min_v = x.floor() as u8;
        let max_v = min_v + 1;
        let frac = x.fract();
        if max_v >= self.data.len() as u8 {
            return *self.data.last().unwrap();
        }
        let up_v: Complex64 = self.data[max_v as usize];
        let down_v: Complex64 = self.data[min_v as usize];
        down_v + ((up_v - down_v) * frac)
    }
}

impl ToFunctionEstimator for &str {
    fn build_estimator<'a>(& self, data: &'a Vec<Complex64>) -> Box<dyn FunctionEstimator<'a> + 'a> {
        match self {
            &"square" => Box::new(SquareInterpolationEstimator::new(data)),
            _ => Box::new(LinearInterpolationEstimator::new(data))
        }
    }
}

pub struct SquareInterpolationEstimator{
    a: Complex64
}

impl <'a> FunctionEstimator<'a> for SquareInterpolationEstimator{
    fn new(data: &'a Vec<Complex64>) -> Self where Self: Sized {
        let real = data.iter().enumerate().map(|(i,s)|i.pow(2) as f64 * (s.re-1.0)).sum::<f64>()
            /
            (0..data.len()).map(|i|i.pow(4) as f64).sum::<f64>();
        let im = data.iter().enumerate().map(|(i,s)|i as f64 * s.im).sum::<f64>()
            /
            (0..data.len()).map(|i|i.pow(2) as f64).sum::<f64>();
        Self{a: Complex64::new(real, im)}
    }

    fn estimate(&self, x: f64) -> Complex64 {
        Complex64::new(1.0 + self.a.re * x.powi(2), self.a.im * x)
    }
}