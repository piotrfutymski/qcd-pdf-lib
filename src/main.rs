use std::path::Path;
use crate::experiment_data::ExperimentData;

mod experiment_data;

fn main() {

    let experiment_data = ExperimentData::load_and_calculate(Path::new("./data/samples"));
    println!("{:?}", experiment_data.reduced_m_value(5, 10).value_error());
}
