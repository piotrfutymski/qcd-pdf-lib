use std::path::Path;
use clap::Parser;
use crate::experiment_data::ExperimentData;

mod experiment_data;

#[derive(Parser)]
struct Args {

    /// Matrix to get data from
    #[arg(short, long)]
    matrix: Option<String>,

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
    num_to_average: Option<u8>,

    /// Type
    #[arg(short, long)]
    polarization_type: Option<String>,

    /// Estimator
    #[arg(short, long)]
    estimator: Option<String>

}

fn main() {

    let args = Args::parse();
    println!("Calculating data from ./data/samples");
    let polarization_type = args.polarization_type.unwrap_or("g0".to_string());
    let mut experiment_data = ExperimentData::load(
        Path::new("./data/samples"),
        args.num_to_average.unwrap_or(8),
        &polarization_type,
        &args.estimator.unwrap_or("linear".to_string())
    );
    println!("Data calculated, printing output");
    let unwrapped_matrix = args.matrix.unwrap_or(String::from("pdf"));
    if unwrapped_matrix.as_str() == "qa"{
        ExperimentData::generate_output_plot_file_one_arg(experiment_data.q_averaged(), args.sample_count.unwrap_or(12), &args.filename.unwrap_or("data.dat".to_string()), !args.imaginary.unwrap_or(false));
    } else if unwrapped_matrix.as_str() == "pdf" {
        let mut experiment_data = ExperimentData::load(
            Path::new("./data/samples"),
            args.num_to_average.unwrap_or(8),
            "g0",
            "linear"
        );
        experiment_data.get_pdf_params_to_file("_g0");
        let mut experiment_data = ExperimentData::load(
            Path::new("./data/samples"),
            args.num_to_average.unwrap_or(8),
            "g0",
            "square"
        );
        experiment_data.get_pdf_params_to_file("_g0_square");

        let mut experiment_data = ExperimentData::load(
            Path::new("./data/samples"),
            args.num_to_average.unwrap_or(8),
            "pol",
            "linear"
        );
        experiment_data.get_pdf_params_to_file("_pol");
        let mut experiment_data = ExperimentData::load(
            Path::new("./data/samples"),
            args.num_to_average.unwrap_or(8),
            "pol",
            "square"
        );
        experiment_data.get_pdf_params_to_file("_pol_square");

        let mut experiment_data = ExperimentData::load(
            Path::new("./data/samples"),
            args.num_to_average.unwrap_or(8),
            "tra",
            "linear"
        );
        experiment_data.get_pdf_params_to_file("_tra");
        let mut experiment_data = ExperimentData::load(
            Path::new("./data/samples"),
            args.num_to_average.unwrap_or(8),
            "tra",
            "square"
        );
        experiment_data.get_pdf_params_to_file("_tra_square");
    } else if unwrapped_matrix.as_str() == "est"{
       experiment_data.generate_estimated_results(args.num_to_average.unwrap())
    } else {
        let data = match unwrapped_matrix.as_str() {
            "m" => experiment_data.reduced_m(),
            "mp" => experiment_data.m_prime(),
            "b" => experiment_data.bare_m(),
            _ => experiment_data.q()
        };
        if let Some(momentum) = args.momentum{
            ExperimentData::generate_output_plot_file(data, momentum, args.sample_count.unwrap_or(12), &args.filename.unwrap_or("data.dat".to_string()), !args.imaginary.unwrap_or(false), unwrapped_matrix == "b");
        }else{
            for i in 0..6 {
                ExperimentData::generate_output_plot_file(data, i, args.sample_count.unwrap_or(12), format!("{}.dat", i).as_str(), !args.imaginary.unwrap_or(false), unwrapped_matrix == "b");
            }
        }

    }

}