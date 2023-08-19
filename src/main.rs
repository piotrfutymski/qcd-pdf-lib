mod bootstrap_data;

use bootstrap_data::BootstrapDataCreator;

fn main() {
    for i in 0..6 {
        BootstrapDataCreator::new(i, "./data/samples", "./data/my_results").create_output().expect("problem with creating output");
    }
}
