mod bootstrap_input;

use std::{fs, io};
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
pub use bootstrap_input::BootstrapInput;

pub struct BootstrapDataCreator{
    mom: i32,
    input_path: String,
    output_path: String
}

impl BootstrapDataCreator {
    pub fn new(mom: i32, input_path: &str, output_path: &str) -> BootstrapDataCreator{
        BootstrapDataCreator{mom, input_path: String::from(input_path), output_path: String::from(output_path)}
    }

    pub fn create_output(&self) -> io::Result<()>{
        let mut inputs: HashMap<i32, BootstrapInput> = self.get_inputs(true)?;
        self.get_inputs(false)?
            .into_iter().for_each(|e|{inputs.insert(e.0, e.1);});
        let mut lines: Vec<(i32, String)> = inputs
            .iter()
            .map(|e|{
                let sample_0 = e.1.get_sample_0();
                let boot_average = e.1.boot_average();
                let boot_error = e.1.boot_error();
                (*e.0, format!("{} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6}\n", e.0, sample_0.re, boot_average.re, boot_error.re, sample_0.im, boot_average.im, boot_error.im))
            })
            .collect();
        lines.sort_by(|a,b|a.0.cmp(&b.0));
        let mut f = File::create(format!("{}/me_mom_{}.dat",self.output_path, self.mom))?;
        lines.iter().for_each(|(_,line)|{
            f.write(line.as_bytes()).unwrap();
        });
        Ok(())
    }

    fn get_inputs(&self, plus_sign: bool) -> io::Result<HashMap<i32, BootstrapInput>> {
        Ok(fs::read_dir(self.input_path.clone())?
            .into_iter()
            .filter_map(|e|e.ok())
            .filter(|e|e.file_name().to_str().unwrap().contains(&format!("mom{}",self.mom)))
            .map(|e|(
                {
                    let num: i32 = e.file_name().to_str().unwrap().split("_").skip(1).next()
                        .unwrap().strip_prefix("z").unwrap().parse().unwrap();
                    match plus_sign {
                        true => num,
                        false => -num
                    }
                }
                ,
                BootstrapInput::read_from_file(e.path().to_str().unwrap(), plus_sign).unwrap())
            )
            .collect())
    }
}
