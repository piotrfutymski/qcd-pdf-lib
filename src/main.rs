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
    imaginary: bool,

    /// Sample count
    #[arg(short, long)]
    sample_count: u8,

    /// Momentum
    #[arg(long)]
    momentum: u8,

    /// Filename
    #[arg(short, long)]
    filename: String

}

fn main() {

    let args = Args::parse();

    println!("Calculating data from ./data/samples");
    let experiment_data = ExperimentData::load_and_calculate(Path::new("./data/samples"));
    println!("Data calculated, printing output");
    let data = match args.matrix.as_str() {
        "m" => experiment_data.reduced_m(),
        "mp" => experiment_data.m_prime(),
        _ => experiment_data.q()
    };
    ExperimentData::generate_output_plot_file(data, args.momentum, args.sample_count, &args.filename, !args.imaginary);

}
