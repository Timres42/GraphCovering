use std::env;
use std::io::Write;
use std::fs::OpenOptions;
use std::time::Instant;
use std::sync::RwLock;
use std::sync::Arc;

pub type Perm = Vec<usize>;
pub type Cycle = Vec<usize>;
pub type CycleSet = Vec<Cycle>;

pub struct CyclePermuter {
    _num_nodes : usize,
    _perm_cycles : Vec<Vec<usize>>,
    _perm_used : Vec<usize>,
    _ham_cycles : Vec<Vec<i16>>,
    _used_edges : Vec<Vec<i8>>,
    _num_ham_cycles : usize,
    _sol_found : Vec<CycleSet>,
}

impl CyclePermuter {

    pub fn new(num_nodes : usize) -> CyclePermuter {
        if(num_nodes-1)%4 != 0 {
            panic!("The number of nodes must be 1 more than a multiple of 4");
        }
        //println!("Creating a CyclePermuter with {} nodes",num_nodes);

        CyclePermuter {
            _num_nodes : num_nodes,
            _perm_cycles : vec![],
            _perm_used : vec![0;num_nodes],
            _ham_cycles : vec![vec![-1;num_nodes];(num_nodes-1)/4],
            _used_edges : vec![vec![0;num_nodes];num_nodes],
            _num_ham_cycles : (num_nodes-1)/4,
            _sol_found : vec![],
        }
    }

    fn inv(&self, i : usize) -> usize{
        return self._num_nodes - i;
    }

    fn set_used(&mut self, i : usize) {
        self._perm_used[i] = 1;
        self._perm_used[self._num_nodes - i] = 1;
    }

    fn set_unused(&mut self, i : usize) {
        self._perm_used[i] = 0;
        self._perm_used[self._num_nodes - i] = 0;
    }

    fn set_ham_cycle_succ(&mut self, cycle_idx : usize, node_idx : usize, val : i16, changed_vals : &mut (Vec<(usize,usize)>,Vec<(usize,usize)>)) -> bool {
        self._ham_cycles[cycle_idx][node_idx] = val;
        changed_vals.0.push((cycle_idx,node_idx));
        //Check for doubly covered edges and update
        let curr_cycl_len = self._ham_cycles[cycle_idx].len();
        let mut edges : Vec<(usize,usize)> = vec![];

        for idx in vec![node_idx+curr_cycl_len-2,node_idx+curr_cycl_len-1,node_idx+1,node_idx+2] {
            let new_val = self._ham_cycles[cycle_idx][idx%curr_cycl_len];
            if new_val != -1 {
                edges.push((new_val as usize,val as usize));
            }
        }
        for edge in edges {
            if self._used_edges[edge.0][edge.1] == 1 {
                return false;
            }
            self._used_edges[edge.0][edge.1] = 1;
            self._used_edges[edge.1][edge.0] = 1;
            changed_vals.1.push((edge.0,edge.1));
        }
        return true;
    }

    fn reset_ham_cycle_from_changed_vals(&mut self, changed_vals : &(Vec<(usize,usize)>,Vec<(usize,usize)>)) {
        for i in 0..changed_vals.0.len() {
            self._ham_cycles[changed_vals.0[i].0][changed_vals.0[i].1] = -1;
        }
        for i in 0..changed_vals.1.len() {
            self._used_edges[changed_vals.1[i].0][changed_vals.1[i].1] = 0;
            self._used_edges[changed_vals.1[i].1][changed_vals.1[i].0] = 0;
        }
    }

    fn reset_ham_cycle(&mut self, cycle_idx : usize, node_idx : usize) {
        if self._ham_cycles[cycle_idx][node_idx] == -1 {
            return;
        }
        let curr_val : usize = self._ham_cycles[cycle_idx][node_idx] as usize;

        self._ham_cycles[cycle_idx][node_idx] = -1;
        //Reset covered edges
        let curr_cycl_len = self._ham_cycles[cycle_idx].len();
        for idx in vec![node_idx+curr_cycl_len-2,node_idx+curr_cycl_len-1,node_idx+1,node_idx+2] {
            let new_val = self._ham_cycles[cycle_idx][idx%curr_cycl_len];
            if new_val != -1 {
                self._used_edges[new_val as usize][curr_val] = 0;
                self._used_edges[curr_val][new_val as usize] = 0;
            }
        }
    }

    fn add_node_to_cycle_succ(&mut self, node : usize, cycle_idx : usize) -> bool{
        //Update perm_cycles
        self._perm_cycles[cycle_idx].push(node);
        //Update perm_used 
        self.set_used(node);

        //Update ham_cycles and used_edges
        let mut change_succ = true;
        let mut changed_values : (Vec<(usize,usize)>,Vec<(usize,usize)>)= (vec![],vec![]);
        for it in 0..self._perm_cycles[cycle_idx].len() -1 {
            let curr_val = self._perm_cycles[cycle_idx][it];
            let curr_cycl_len = self._perm_cycles[cycle_idx].len();
            change_succ = 
            self.set_ham_cycle_succ(curr_cycl_len - 1 - it, curr_val, node as i16, &mut changed_values) 
            && self.set_ham_cycle_succ(self._num_ham_cycles + 1 + it - curr_cycl_len, node, self.inv(curr_val) as i16, &mut changed_values)
            && self.set_ham_cycle_succ(curr_cycl_len - 1 - it, self.inv(curr_val), self.inv(node) as i16, &mut changed_values)
            && self.set_ham_cycle_succ(self._num_ham_cycles + 1 + it - curr_cycl_len, self.inv(node), curr_val as i16, &mut changed_values);

            if !change_succ {
                break;
            }
        }
        if !change_succ {
            self.reset_ham_cycle_from_changed_vals(&changed_values);
            self.set_unused(node);
            self._perm_cycles[cycle_idx].pop();
            return false;
        }
        return true;
    }

    fn remove_node_from_cycle(&mut self, node : usize, cycle_idx : usize) {
        //Update ham_cycles and used_edges
        let curr_cycl_len = self._perm_cycles[cycle_idx].len();
        for it in 0..self._perm_cycles[cycle_idx].len() -1 {
            let curr_val = self._perm_cycles[cycle_idx][it];
            self.reset_ham_cycle(curr_cycl_len - 1 - it, curr_val);
            self.reset_ham_cycle(self._num_ham_cycles + 1 + it - curr_cycl_len, node);
            self.reset_ham_cycle(curr_cycl_len - 1 - it, self.inv(curr_val));
            self.reset_ham_cycle(self._num_ham_cycles + 1 + it - curr_cycl_len, self.inv(node));
        }
        //Update perm_cycles
        self._perm_cycles[cycle_idx].pop();
        //Update perm_used
        self.set_unused(node);
    }

    pub fn enum_all_symm_perms(&mut self, first_node : usize) {
        self._perm_cycles = vec![vec![1],vec![2]];
        self.set_used(1);
        self.set_used(2);
        if self.add_node_to_cycle_succ(first_node, 0) {
            self.rec_enum_symm_perms();
            self.remove_node_from_cycle(first_node, 0);
        }
    }

    fn rec_enum_symm_perms(&mut self) {
        if self._perm_cycles[1].len() == (self._num_nodes-1)/4 {
            //Add curr perm and cycles to sol_found
            self._sol_found.push(self._perm_cycles.clone());
            return;
        }
        
        //let upper_bound = if self._perm_cycles[0].len() == 1 {(self._num_nodes - 1)/2} else {self._num_nodes -2};
        let upper_bound = self._num_nodes - 2;
        for node in 3..upper_bound {
            if self._perm_used[node] == 1 {continue;}
            let cycle_idx = if self._perm_cycles[0].len() < self._num_ham_cycles {0} else {1};
            if self.add_node_to_cycle_succ(node, cycle_idx) {
                self.rec_enum_symm_perms();
                self.remove_node_from_cycle(node, cycle_idx);
            }
        }
    }
}

fn run_threaded_cycle_permuter(num_nodes : usize, num_threads : usize) {
    println!("Running CyclePermuter with {} nodes and {} threads",num_nodes,num_threads);
    //Initialize thread pool
    let pool = rayon::ThreadPoolBuilder::new().num_threads(num_threads).build().unwrap();
    let lock = Arc::new(RwLock::new(vec![]));
    let lock_owned = lock.clone();

    let start = Instant::now();
    //Start multithreaded computation
    pool.scope(move |s| {
        //Create a thread for each permutation
        for first_node in 3..num_nodes-2 {
            let lock_owned = lock_owned.clone();
            s.spawn(move |_| {
                let mut cp = CyclePermuter::new(num_nodes);
                cp.enum_all_symm_perms(first_node);
                //println!("Thread {} finished with {} number of solutions.",first_node,cp._sol_found.len());

                let mut result = lock_owned.write().unwrap();
                for sol in cp._sol_found {
                    result.push(sol);
                }
                drop(result);
            });
        }
    });

    let result = lock.read().unwrap();
    println!("Number of solutions found: {}",result.len());
    let filename = "sol/perms_".to_owned() + &num_nodes.to_string() + ".txt";
        let mut _file = OpenOptions::new().write(true).truncate(true).create(true).open(filename).unwrap();
    for sol in result.iter() {
        let mut perm_cycles_copy = sol.clone();
        //Add inverses to each cycle to get the full permutation cycles
        for cycle in &mut perm_cycles_copy {
            let inverses: Vec<_> = cycle.iter().map(|&node| num_nodes - node).collect();
            cycle.extend(inverses);
        }
        //Write cycles to file
        for cycle in &perm_cycles_copy {
            let _ = write!(_file, "({}) ", cycle.iter().map(|node| node.to_string()).collect::<Vec<_>>().join(" "));
        }
        let _ = writeln!(_file);
    }
    drop(result);
    let duration = start.elapsed();
    println!("Time elapsed is: {:?}", duration);
    println!();
}

//Main has one argument
fn main () {
    //Read in command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        println!("Usage: {} <num_threads>", args[0]);
        return;
    }
    let num_threads = args[1].parse::<usize>().unwrap();

    for num_nodes in vec![21,25,29,33,37] {
        run_threaded_cycle_permuter(num_nodes,num_threads);
    }    
}
