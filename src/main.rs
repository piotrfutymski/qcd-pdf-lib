use std::path::Path;
use clap::Parser;
use crate::experiment_data::ExperimentData;

mod experiment_data;

#[derive(Parser)]
struct Args {

    /// Matrix to get data from
    #[arg(short, long)]
    matrix: String,

    /// Data from imaginary Part
    #[arg(short, long)]
    imaginary: Option<bool>,

    /// Sample count
    #[arg(short, long)]
    sample_count: Option<u8>,

    /// Momentum
    #[arg(long)]
    momentum: Option<u8>,

    /// Filename
    #[arg(short, long)]
    filename: Option<String>,

    /// Num to average
    #[arg(short, long)]
    num_to_average: Option<u8>

}

fn main() {

    let args = Args::parse();

    println!("Calculating data from ./data/samples");
    let mut experiment_data = ExperimentData::load(Path::new("./data/samples"), args.num_to_average.unwrap_or(8));
    println!("Data calculated, printing output");
    if args.matrix.as_str() == "qa"{
        ExperimentData::generate_output_plot_file_one_arg(experiment_data.q_averaged(), args.sample_count.unwrap_or(23), &args.filename.unwrap_or("data.dat".to_string()), !args.imaginary.unwrap_or(false));
    } else {
        let data = match args.matrix.as_str() {
            "m" => experiment_data.reduced_m(),
            "mp" => experiment_data.m_prime(),
            _ => experiment_data.q()
        };
        ExperimentData::generate_output_plot_file(data, args.momentum.unwrap_or(0), args.sample_count.unwrap_or(23), &args.filename.unwrap_or("data.dat".to_string()), !args.imaginary.unwrap_or(false));
    }

}


#[test]
fn test() {
    let experiment_data = ExperimentData::load(Path::new("./data/samples"), 8);
}