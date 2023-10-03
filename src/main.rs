mod bootstrap_data_creator;
mod bootstrap_data;
mod experiment_data;

use bootstrap_data_creator::BootstrapDataCreator;

fn main() {
    for i in 0..6 {
        BootstrapDataCreator::new(i, "./data/samples", "./data/my_results").create_output().expect("problem with creating output");
    }
}
