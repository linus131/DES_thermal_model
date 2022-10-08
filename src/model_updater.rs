use crate::model_generator::ModelGenerator;
use std::collections::HashMap;

/// ModelUpdater stores information about which cell is activated at what time
#[derive(Clone)]
pub struct ModelUpdater{
    pub(crate) activation_times: Vec<(f64, usize, [f64;2], usize)>,
    pub(crate) model_generator: ModelGenerator,
    is_cell_activated: Vec<bool>,
    current_cell:usize,
    nodenum_new: HashMap<usize,usize>
}
/// methods for ModelUpdater
impl ModelUpdater{
    pub fn new(activation_times: Vec<(f64,usize,[f64;2],usize)>, model_generator: ModelGenerator, nodenum_new: HashMap<usize,usize>)->ModelUpdater{
        let is_cell_activated = (0..model_generator.celllist.len()).map(|i| false).collect();
        return ModelUpdater{
            activation_times,
            model_generator,
            is_cell_activated,
            current_cell:0,
            nodenum_new,
        }
    }

    pub fn get_next_cell_info(&mut self)->(f64,[usize;8],[isize;6],[bool;6],usize,[f64;2],usize){
        let cell_no = self.activation_times[self.current_cell].1;
        let activation_time = self.activation_times[self.current_cell].0;
        let orientation = self.activation_times[self.current_cell].2;
        let layer_no = self.activation_times[self.current_cell].3;


        let neighborinfo = self.model_generator.neighborlist[cell_no];
        let mut neighbor_bool = self.model_generator.has_neighbor[cell_no];
        self.is_cell_activated[cell_no] = true;
        if neighbor_bool[0]{
            if !self.is_cell_activated[neighborinfo[0] as usize] { neighbor_bool[0] = false}
        }
        if neighbor_bool[1]{
            if !self.is_cell_activated[neighborinfo[1] as usize] { neighbor_bool[1] = false}
        }
        if neighbor_bool[2]{
            if !self.is_cell_activated[neighborinfo[2] as usize] { neighbor_bool[2] = false}
        }
        if neighbor_bool[3]{
            if !self.is_cell_activated[neighborinfo[3] as usize] { neighbor_bool[3] = false}
        }
        if neighbor_bool[4]{
            if !self.is_cell_activated[neighborinfo[4] as usize] { neighbor_bool[4] = false}
        }
        if neighbor_bool[5]{
            if !self.is_cell_activated[neighborinfo[5] as usize] { neighbor_bool[5] = false}
        }
        self.current_cell = self.current_cell + 1;

        let mut newelem = [0;8];
        for i in 0..8{
            newelem[i] = *self.nodenum_new.get(&self.model_generator.celllist[cell_no][i]).expect("can't find the old node to new node mapping");
        }
        return(activation_time,newelem,neighborinfo,neighbor_bool,cell_no,orientation,layer_no);

    }


}
