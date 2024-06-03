use aseq_lib::{distance_map, filter_seq_greedy, filter_seq_hybrid, get_all_sequences};
use std::{
    env,
    io::{Read, Write},
};

#[derive(Debug)]
enum Filter {
    Greedy,
    Hybrid(f64),
}

fn main() {
    // Load vectors from given file
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        println!(
            r#"
Invalid number of arguments.
Usage: aseq <path_to_vectors_file> <filter> [options]
Where filter can be:
  'greedy' - no additional options
  'hybrid' - [cutoff] - value between 0 and 1. Sequences with coverage to lenght ratio
                        below cutoff will be discarded during filtering
"#
        );
        return;
    }

    let filter = match args[2].as_str() {
        "greedy" => Filter::Greedy,
        "hybrid" => {
            if args.len() >= 4 {
                let cutoff: f64 = args[3].parse().expect("");
                Filter::Hybrid(cutoff)
            } else {
                Filter::Hybrid(0.5)
            }
        }
        _ => {
            println!("Invalid filter. Valid filters are: 'greedy', 'hybrid'");
            return;
        }
    };

    let file = std::fs::File::open(&args[1]);
    // let file = std::fs::File::open("test.txt");
    let mut content = String::new();
    let _ = match file {
        Ok(mut f) => f.read_to_string(&mut content),
        Err(e) => {
            println!("Error opening file: {}", e);
            return;
        }
    };

    let mut vectors: Vec<u64> = Vec::new();
    for line in content.lines() {
        let v: u64 = line
            .parse()
            .expect(&format!("Couldn't parse line: {} to u64", line));
        if vectors.contains(&v) {
            println!("Vectors must be unique. Duplicate item: {}", v);
            return;
        }
        vectors.push(v);
    }

    // create distance map
    let map = distance_map(&vectors);

    let mut csv = String::new();
    for row in map.iter() {
        for item in row.iter() {
            csv.push_str(&format!("{},", item));
        }
        csv.push('\n');
    }
    let mut csv_file = std::fs::File::create("distances.csv").unwrap();
    let _ = csv_file.write_all(&csv.as_bytes());
    println!("Distance map saved to distances.csv file.");

    // create sequences
    let res = get_all_sequences(&map);
    println!("Sequences");
    println!("{:?}", res);
    println!("Sequences de-mapped");
    for seq in res.iter() {
        for idx in seq[0..seq.len() - 1].iter() {
            print!("{} -> ", vectors[*idx]);
        }
        println!("{}", vectors[seq[seq.len() - 1]]);
    }

    // Create pattern set
    let patterns = match filter {
        Filter::Hybrid(cutoff) => filter_seq_hybrid(res, &vectors, cutoff),
        Filter::Greedy => filter_seq_greedy(res, &vectors),
    };
    println!("Filtered patterns: {:?}", patterns);
}
