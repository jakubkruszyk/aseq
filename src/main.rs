use std::{
    collections::HashSet,
    env,
    io::{Read, Write},
};

fn main() {
    // Load vectors from given file
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        println!("Invalid number of arguments.\nUsage: aseq <path_to_vectors_file>");
        return;
    }

    let file = std::fs::File::open(&args[1]);
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
        csv.push_str("\n");
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
    let patterns = filter_seq_bms(res, &vectors);
    println!("Filtered patterns: {:?}", patterns);
}

/// Creates distance map
fn distance_map(vectors: &[u64]) -> Vec<Vec<u64>> {
    let mut map: Vec<Vec<u64>> = Vec::with_capacity(vectors.len());
    for item in vectors.iter() {
        let mut row: Vec<u64> = Vec::with_capacity(vectors.len());
        for next_item in vectors.iter() {
            let diff = *next_item as i128 - *item as i128;
            let udiff = if diff < 0 {
                (!diff + 1) as u64
            } else {
                diff as u64
            };
            row.push(udiff);
        }
        map.push(row);
    }
    map
}

/// Create arithmetic sequence based on given element (row) and next element (col)
fn get_single_seq(map: &Vec<Vec<u64>>, row: usize, col: usize) -> Vec<usize> {
    // vector of number ids that forms arithmetic sequence
    let mut path: Vec<usize> = vec![row];
    let increment = map[row][col];
    let mut next = col;
    while !path.contains(&next) && path.len() < map.len() {
        path.push(next);
        for (idx, item) in map[next].iter().enumerate() {
            if *item == increment {
                next = idx;
            }
        }
    }
    path
}

/// Check if seq2 is part of seq1
fn in_seq(seq1: &[usize], seq2: &[usize]) -> bool {
    if seq2.len() > seq1.len() {
        return false;
    }
    let first_match = seq1.iter().position(|x| *x == seq2[0]);
    match first_match {
        Some(pos) => {
            // check if there is enough space in seq1 slice to contain the seq2
            if seq2.len() > seq1.len() - pos {
                return false;
            }
            // compare sequences
            let mut res = true;
            for (s1, s2) in seq2.iter().zip(seq1[pos..].iter()) {
                if *s1 != *s2 {
                    res = false;
                }
            }
            res
        }
        None => false,
    }
}

/// Generate all posible sequences with length >= 3 (no duplicates)
fn get_all_sequences(map: &Vec<Vec<u64>>) -> Vec<Vec<usize>> {
    let mut sequences: Vec<Vec<usize>> = Vec::new();
    for (r_idx, row) in map.iter().enumerate() {
        for (c_idx, item) in row.iter().enumerate() {
            // skip when next row is the same as current row
            if *item == 0 {
                continue;
            }
            let seq = get_single_seq(map, r_idx, c_idx);
            if seq.len() >= 3 && !sequences.contains(&seq) {
                sequences.push(seq);
            }
        }
    }
    sequences
}

/// Remove all items present in 'seq' from other sequences passed as 'sequences'
/// If 'divide' is false, and item was found inside a sequence, it will be discarded
/// otherwise it will be split in half and only item will be removed
fn eliminate_sequence(
    mut sequences: Vec<Vec<usize>>,
    seq: &Vec<usize>,
    divide: bool,
) -> Vec<Vec<usize>> {
    let mut filtered_seq: Vec<Vec<usize>> = Vec::new();
    // treat 'sequences' vector as queue of sequences that need to be checked
    'outer: while sequences.len() > 0 {
        let mut s = sequences.pop().unwrap();
        for item in seq.iter() {
            if let Some(pos) = s.iter().position(|x| *x == *item) {
                if pos == 0 || pos == s.len() - 1 {
                    s.remove(pos);
                    if s.len() > 0 {
                        sequences.push(s);
                    }
                } else if divide {
                    let s2 = s.split_off(pos + 1);
                    let _ = s.pop();
                    if s.len() > 0 {
                        sequences.push(s);
                    }
                    if s2.len() > 0 {
                        sequences.push(s2);
                    }
                }
                continue 'outer;
            };
        }
        // Filter sequences shorter than 3 and duplicates
        if s.len() > 2 && !filtered_seq.contains(&s) {
            filtered_seq.push(s);
        }
    }
    filtered_seq
}

/// Filter sequences to cover all vectors using sequences that matches the most of remaining
/// symbols
fn filter_seq_bms(mut sequences: Vec<Vec<usize>>, vectors: &Vec<u64>) -> Vec<Vec<usize>> {
    let mut symbols: HashSet<usize> = HashSet::new();
    for i in 0..vectors.len() {
        symbols.insert(i);
    }
    let mut filtered_seq: Vec<Vec<usize>> = Vec::new();
    while symbols.len() > 0 && sequences.len() > 0 {
        // Get sequence that has the highest ratio of covered remaining symbols per sequence length
        let mut idx = 0;
        let mut ratio = 0.0;
        for (i, seq) in sequences.iter().enumerate() {
            let r = items_in_seq(&seq, &symbols) as f64 / seq.len() as f64;
            if r > ratio {
                ratio = r;
                idx = i;
            }
        }
        let seq = sequences.remove(idx);
        for item in seq.iter() {
            symbols.remove(item);
        }
        filtered_seq.push(seq);
    }
    // Add remaining sequences as single vectors
    for symbol in symbols.into_iter() {
        filtered_seq.push(vec![symbol]);
    }
    filtered_seq
}

fn items_in_seq(seq: &Vec<usize>, symbols: &HashSet<usize>) -> usize {
    let mut count = 0;
    for item in seq.iter() {
        if symbols.contains(item) {
            count += 1;
        }
    }
    count
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_distance_map() {
        // let test_seq = [1, 2, 5, 4, 3, 6];
        // let map = distance_map(&test_seq);
        // let map_gold = [
        //     [0, 1, 4, 3, 2, 5],
        //     [-1, 0, 3, 2, 1, 4],
        //     [-4, -3, 0, -1, -2, 1],
        //     [-3, -2, 1, 0, -1, 2],
        //     [-2, -1, 2, 1, 0, 3],
        //     [-5, -4, -1, -2, -3, 0],
        // ];
        // for (row, row_g) in map.iter().zip(map_gold.iter()) {
        //     for (item, item_g) in row.iter().zip(row_g.iter()) {
        //         assert_eq!(*item, *item_g);
        //     }
        // }
    }

    #[test]
    fn test_get_sequence() {
        let map = distance_map(&[1, 2, 5, 4, 3, 6]);
        let seq = get_single_seq(&map, 0, 1);
        let seq_g = [0, 1, 4, 3, 2, 5];
        assert_eq!(seq.len(), seq_g.len());
        for (item, item_g) in seq.iter().zip(seq_g) {
            assert_eq!(*item, item_g);
        }
        let seq = get_single_seq(&map, 0, 4);
        let seq_g = [0, 4, 2];
        assert_eq!(seq.len(), seq_g.len());
        for (item, item_g) in seq.iter().zip(seq_g) {
            assert_eq!(*item, item_g);
        }
    }

    #[test]
    fn test_in_seq() {
        let seq1 = [1, 2, 3, 4, 5, 6];
        let seq2 = [1, 2, 3];
        assert_eq!(in_seq(&seq1, &seq2), true);
        let seq2 = [3, 4, 5];
        assert_eq!(in_seq(&seq1, &seq2), true);
        let seq2 = [5, 6, 7];
        assert_eq!(in_seq(&seq1, &seq2), false);
        let seq2 = [1, 2, 4, 5];
        assert_eq!(in_seq(&seq1, &seq2), false);
        let seq2 = [1, 2, 3, 4, 5, 6, 7];
        assert_eq!(in_seq(&seq1, &seq2), false);
    }
}
