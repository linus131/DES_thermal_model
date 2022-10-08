use crate::primitives::Node;
use std::collections::HashMap;
use crate::interpolator::Interpolator;
use rayon::ThreadPool;
use std::io::BufWriter;
use std::fs::File;
use crate::model_updater::ModelUpdater;
use std::io::Write;

/// Model has list of Nodes, list of cells, list of temperature for each cell, and information about
/// six neighbors (other cell or free surface) for each cell
pub struct Model{
    nodelist:Vec<Node>,
    sp_heat_cap:Vec<f64>,
    density:Vec<f64>,
    kx:Vec<f64>,
    pub(crate) templist:Vec<f64>,
    has_neighbor:Vec<[bool;6]>,
    neighbors:Vec<[isize;6]>,
    pub(crate) temp_at_contact:Vec<[f64;6]>,
    density_vol_sp_ht_cap:Vec<f64>,
    pub(crate) cell_index_map:HashMap<isize,isize>,
    sp_heat_cap_interp:Interpolator,
    element_width:f64,
    element_height:f64,
    active_layer_first_element:usize,
    previous_active_layer:usize,
    previous_layer:usize,
    current_layer:usize,
    turn_off_layers_at:usize,
    first_element_of_each_layer:Vec<usize>,
    temp_pos_for_split:Vec<usize>,
}

/// methods for Model struct
impl Model{
    /// generates a new empty Model
    pub fn new(nodes:Vec<Node>,activation_cells_length:usize, sp_heat_cap_interp: Interpolator, element_width:f64, element_height:f64,
               turn_off_layers_at:usize)->Model{
        return Model{
            nodelist: nodes,
            sp_heat_cap:Vec::with_capacity(activation_cells_length),
            density:Vec::with_capacity(activation_cells_length),
            kx:Vec::with_capacity(activation_cells_length),
            templist: Vec::with_capacity(activation_cells_length),
            neighbors: Vec::with_capacity(activation_cells_length),
            has_neighbor:Vec::with_capacity(activation_cells_length),
            temp_at_contact: Vec::with_capacity(activation_cells_length),
            density_vol_sp_ht_cap:Vec::with_capacity(activation_cells_length),
            cell_index_map:HashMap::with_capacity(activation_cells_length),
            sp_heat_cap_interp,
            element_width,
            element_height,
            active_layer_first_element:0,
            previous_active_layer:0,
            previous_layer: 0,
            current_layer: 0,
            turn_off_layers_at,
            first_element_of_each_layer: vec![0],
            temp_pos_for_split:Vec::with_capacity(activation_cells_length),
        }
    }
    /// adds a cell to the model. Takes node index, specific heat capacity, density, conductivity in
    /// the x-direction, conductivity in the y and z-directions, initial temperature, and neighbor
    /// information of the cell
    pub fn addCell(&mut self,node_no:[usize;8],sp_heat_cap:f64,
                   density:f64, kx:f64, ky:f64, kz:f64, init_temp:f64,is_neighbor:[bool;6],
                   neighbor:[isize;6],cell_index:usize,temporary_templist:&mut Vec<f64>,
                   orientation:[f64;2],layer:usize,cell_indices:&mut Vec<usize>)
    {
        let mut counter = self.sp_heat_cap.len();
        self.cell_index_map.insert(cell_index as isize,counter as isize);
        cell_indices.push(counter);
        if (layer != self.current_layer){
            self.previous_layer = self.current_layer;
            self.current_layer = layer;
            self.first_element_of_each_layer.push(counter);
            //println!("layer {}",layer);
        }
        if (layer-self.previous_active_layer)>self.turn_off_layers_at{
            self.previous_active_layer = self.previous_active_layer+1;
            self.active_layer_first_element = self.first_element_of_each_layer[self.previous_active_layer+1];

        }

        let p1 = self.nodelist[node_no[0]].clone();
        let p2 = self.nodelist[node_no[7]].clone();

        let volume = self.element_width*self.element_width*self.element_height;//((p2.pt.x-p1.pt.x)*(p2.pt.y-p1.pt.y)*(p2.pt.z-p1.pt.z)).abs();

        self.density_vol_sp_ht_cap.push(density*volume);
        temporary_templist.push(init_temp.clone());


        // please check this one more time
        self.kx.push(kx);
        // self.ky.push(ky);

        //self.kz.push(kz);
        self.density.push(density);
        self.sp_heat_cap.push(sp_heat_cap);
        // println!("start suspicious code");

        //update neighbor index
        let mut neighbor2:[isize;6] =[0,0,0,0,0,0];
        self.temp_at_contact.push([-100.0;6]);
        for i in 0..6{
            if is_neighbor[i]{
                neighbor2[i] = self.cell_index_map.get(&(neighbor[i])).expect("can't find it").clone();
            }
        }

        self.templist.push(init_temp);
        if is_neighbor[0] {
            self.neighbors[neighbor2[0] as usize][1] = counter as isize;
            self.has_neighbor[neighbor2[0] as usize][1] = true;
            self.temp_at_contact[neighbor2[0] as usize][1] = self.templist[neighbor2[0] as usize];
            //self.temp_at_contact[neighbor2[0] as usize][1] = self.templist[neighbor2[0] as usize][1]
            //  self.has_neighbor_f64[neighbor2[0] as usize][1] = 1.0;
            //  self.has_no_neighbor_f64[neighbor2[0] as usize][1] = 0.0;
        }
        //println!("check 1 pass");
        if is_neighbor[1] {
            self.neighbors[neighbor2[1] as usize][0]=counter as isize;
            self.temp_at_contact[neighbor2[1] as usize][0] = self.templist[neighbor2[1] as usize];
            self.has_neighbor[neighbor2[1] as usize][0]=true;
            // self.has_neighbor_f64[neighbor2[1] as usize][0] = 1.0;
            //self.has_no_neighbor_f64[neighbor2[1] as usize][0] = 0.0;
        }
        //  println!("check 2 pass");
        if is_neighbor[2] {
            self.neighbors[neighbor2[2] as usize][3]=counter as isize;
            self.temp_at_contact[neighbor2[2] as usize][3] = self.templist[neighbor2[2] as usize];
            self.has_neighbor[neighbor2[2] as usize][3]=true;
            // self.has_neighbor_f64[neighbor2[2] as usize][3] = 1.0;
            // self.has_neighbor_f64[neighbor2[2] as usize][3] = 0.0;
        }
        // println!("check 3 pass");
        if is_neighbor[3] {
            self.neighbors[neighbor2[3] as usize][2]=counter as isize;
            self.temp_at_contact[neighbor2[3] as usize][2] = self.templist[neighbor2[3] as usize];
            self.has_neighbor[neighbor2[3] as usize][2] = true;
            //self.has_neighbor_f64[neighbor2[3] as usize][2] = 1.0;
            //self.has_neighbor_f64[neighbor2[3] as usize][2] = 0.0;
        }
        // println!("check 4 pass");
        if is_neighbor[4] {
            self.neighbors[neighbor2[4] as usize][5]=counter as isize;
            self.temp_at_contact[neighbor2[4] as usize][5] = self.templist[neighbor2[4] as usize];
            self.has_neighbor[neighbor2[4] as usize][5]=true;
        }

        if is_neighbor[5] {
            self.neighbors[neighbor2[5] as usize][4]= counter as isize;
            self.temp_at_contact[neighbor2[5] as usize][4] = self.templist[neighbor2[5] as usize];
            self.has_neighbor[neighbor2[5] as usize][4] = true;
        }

        self.has_neighbor.push(is_neighbor);
        self.neighbors.push(neighbor2);

    }


    /// calculates new cell temperatures for the model for a timestep. Takes time increment,
    /// convection coefficient and environment temperature as input. Gives out temperature of the
    /// cells as output

    pub fn find_new_cell_temp_all(&self, dt:f64, conv_coeff:f64, t_env:f64, side_area:f64, top_area:f64, sa_d_sd:f64, ta_d_td:f64,
                                  vol:f64,pool:&ThreadPool,
                                  newtemp:&mut [f64],numthreads:usize, maxthreads:usize, cell_indices_slice: &[usize]) {
        if (numthreads < maxthreads){
            let newtempsize = newtemp.len();
            let (nt0,nt1) = newtemp.split_at_mut(newtempsize/2);
            let (cellpos0, cellpos1) = cell_indices_slice.split_at(newtempsize/2);
            pool.install(
                ||rayon::join(
                    || self.find_new_cell_temp_all(dt, conv_coeff,t_env,side_area,top_area, sa_d_sd,ta_d_td,vol, pool, nt0,numthreads*2,maxthreads,cellpos0),
                    || self.find_new_cell_temp_all(dt, conv_coeff,t_env,side_area,top_area, sa_d_sd,ta_d_td,vol, pool, nt1,numthreads*2,maxthreads,cellpos1)
                ));

        }
        else{
            self.find_new_temp_all_split(dt, conv_coeff, t_env, newtemp, side_area, top_area, sa_d_sd, ta_d_td, vol, cell_indices_slice);
        }
    }



    ///finds cell temperature for values in a range
    pub fn find_new_temp_all_split(&self, dt:f64, conv_coeff:f64, t_env:f64, temp: &mut [f64], side_area:f64, top_area:f64, sa_d_sd:f64, ta_d_td:f64, vol:f64,
                                   cell_indices_slice:&[usize]){
        for cell in cell_indices_slice{
            let mut dQ = 0.0;
            //front back
            //dist area and vol is  [ side area, top area, sidearea/sidedistance, toparea/topdistance, volume]
            for i in 0..2 {
                if !self.has_neighbor[*cell][i] {
                    dQ = dQ - conv_coeff * (self.templist[*cell] - t_env) * dt * side_area;// * self.has_no_neighbor_f64[*cell][i] ;
                } else {
                    let neighbor_index = self.neighbors[*cell][i];
                    let t1 = self.templist[*cell];
                    let t2 = self.templist[neighbor_index as usize];
                    dQ = dQ - self.kx[*cell] * (t1 - t2) * sa_d_sd * dt;// *self.has_neighbor_f64[*cell][i];
                }
            }
            //right left
            for i in 2..4 {
                if !self.has_neighbor[*cell][i] {
                    dQ = dQ - conv_coeff * (self.templist[*cell] - t_env) * dt * side_area;// * self.has_no_neighbor_f64[*cell][i] ;
                } else {
                    let neighbor_index = self.neighbors[*cell][i];
                    let t1 = self.templist[*cell];
                    let t2 = self.templist[neighbor_index as usize];
                    dQ = dQ - self.kx[*cell] * (t1 - t2) * sa_d_sd * dt;// *self.has_neighbor_f64[*cell][i];
                }
            }

            //top down
            for i in 4..6 {
                if !self.has_neighbor[*cell][i] {
                    dQ = dQ - conv_coeff * (self.templist[*cell] - t_env) * dt * top_area;
                } else {
                    let neighbor_index = self.neighbors[*cell][i];
                    let t1 = self.templist[*cell];
                    let t2 = self.templist[neighbor_index as usize];
                    dQ = dQ - self.kx[*cell] * (t1 - t2) * ta_d_td * dt;
                }
            }
            let mult = &self.density_vol_sp_ht_cap[*cell] * self.sp_heat_cap[*cell];


            temp[*cell-cell_indices_slice[0]] = (mult * self.templist[*cell] + dQ)/ mult;
            //let newtemperature = temp[*cell-cell_indices_slice[0]];
            /*if newtemperature>220.0 || newtemperature<25.0{
                println!("newtemperature {} *cell-cell_indices_slice[0] {} mult {}", newtemperature, *cell-cell_indices_slice[0], mult)
            }*/
        }
    }


    /// calculates the temperature of the cells in the model during a time period with a given time
    /// step. Takes time period, time step, and a file (buffered) where output is to be written as
    /// input. Writes output to the file (buffer). Might need to change that for speed.
    pub fn run_model(&mut self, time:f64, dt:f64, conv_coeff:f64, t_env:f64, file:&mut BufWriter<File>,
                     global_time:&mut f64,pool:&ThreadPool, areas_and_dists:[f64;5], newtemp: &mut Vec<f64>,
                     conductivity_interp:&Interpolator, maxthreads:usize, cell_indices:&Vec<usize>) {

        let count = (time / dt) as usize;

        let (inactive_temps,active_temps) = newtemp.split_at_mut(self.active_layer_first_element);
        let (old_cell_indices_slice, new_cell_indices_slice) = cell_indices.split_at(self.active_layer_first_element);
        let new_cell_slice_start = new_cell_indices_slice[0].clone();
        let ncis_len = new_cell_indices_slice.len();
        let new_cell_slice_end = new_cell_indices_slice[ncis_len-1].clone();
        {
            let (inactive_selftemps, active_selftemps)  = self.templist.split_at(self.active_layer_first_element);
            let (sp_ht_cap_inactive, sp_ht_cap_active) = self.sp_heat_cap.split_at_mut(self.active_layer_first_element);
            self.sp_heat_cap_interp.interpolate_slice(active_selftemps, sp_ht_cap_active, pool);

            let (inactive_conductivity, active_conductivity) = self.kx.split_at_mut(self.active_layer_first_element);
            conductivity_interp.interpolate_slice(active_selftemps,active_conductivity, pool);

        }



        for i in 0..count {
            self.find_new_cell_temp_all(dt, conv_coeff, t_env, areas_and_dists[0], areas_and_dists[1], areas_and_dists[2],
                                        areas_and_dists[3], areas_and_dists[4], &pool,
                                        active_temps,1,maxthreads, new_cell_indices_slice);

            self.templist[new_cell_slice_start..new_cell_slice_end+1].copy_from_slice(active_temps);
            //if (self.templist.len()>0) {write!(file, "{},{},{},{}\n", self.templist[0].clone(), global_time,self.sp_heat_cap[0],self.density_vol_sp_ht_cap[0] );}

        }
        self.find_new_cell_temp_all(time - count as f64 * dt, conv_coeff, t_env, areas_and_dists[0], areas_and_dists[1], areas_and_dists[2],
                                    areas_and_dists[3], areas_and_dists[4], &pool,
                                    active_temps,1,maxthreads, new_cell_indices_slice);
        self.templist[new_cell_slice_start..new_cell_slice_end+1].copy_from_slice(active_temps);


        if (newtemp.len()>62) {write!(file, "{},{}\n", newtemp[62].clone(), global_time); }


        //12886
        //12473
    }


    /// updates model based on input. Needs to be placed in Model implementation later.
    pub fn update_model(mu:&mut ModelUpdater,mdl:&mut Model,bw:&mut BufWriter<File>,areas_and_dists:[f64;5],time_step:f64,
                        conductivity_interp : &Interpolator,maxthreads:usize,input_data:([f64;3],f64,f64,f64,[f64;2]),pool:&ThreadPool){
        //let max_cpus = num_cpus::get_physical();

        //let pool = rayon::ThreadPoolBuilder::new().num_threads(maxthreads).build().unwrap();
        let mut temporary_templist = Vec::with_capacity(mu.activation_times.len());
        let mut cell_indices = Vec::with_capacity(mu.activation_times.len());
        //println!("activation times len {}" ,mu.activation_times.len());
        let mut next_cell_info = mu.get_next_cell_info();
        let mut next_cell_info2= mu.get_next_cell_info();
        let mut global_time = 0.0;
        let time = next_cell_info2.0 - next_cell_info.0;
        let dt = time_step;
        let kx = input_data.0[0];
        let ky = input_data.0[1];
        let kz = input_data.0[2];
        let density = input_data.1;
        let init_temp = input_data.4[0];
        let t_env = input_data.4[1];
        let conv_coeff = input_data.3;
        let sp_heat_cap = input_data.2;
        mdl.addCell(next_cell_info.1,sp_heat_cap,density,kx,ky,kz,init_temp,
                    next_cell_info.3,next_cell_info.2,next_cell_info.4, &mut temporary_templist,
                    next_cell_info.5,next_cell_info.6, &mut cell_indices);
        //global_time = global_time + time;
        global_time = global_time+time;
        mdl.run_model(time,dt,conv_coeff,t_env,bw, &mut global_time,&pool,areas_and_dists, &mut temporary_templist,
                      conductivity_interp,maxthreads,&cell_indices);

        //println!("run model");
        let endval = mu.activation_times.len();
        for i in 1..endval-1 {

            next_cell_info = next_cell_info2;
            next_cell_info2 = mu.get_next_cell_info();

            let time = next_cell_info2.0 - next_cell_info.0;
            // global_time = global_time + time;
            if i%1000 == 0 {println!("{} out of {}, time {}",i, endval, global_time);}
            //println!("layer {}",next_cell_info.6);
            mdl.addCell(next_cell_info.1,sp_heat_cap,density,kx,ky,kz,init_temp,
                        next_cell_info.3,next_cell_info.2,next_cell_info.4, &mut temporary_templist,
                        next_cell_info.5,next_cell_info.6, &mut cell_indices);
            //println!("second addcell");
            //println!("second run_model");
            global_time = global_time+time;
            if time > 1e-4 {
                mdl.run_model(time, dt, conv_coeff, t_env, bw, &mut global_time, &pool, areas_and_dists, &mut temporary_templist,
                              conductivity_interp, maxthreads, &cell_indices);
            }

            //println!("second run model");
        }
        //println!("pulupulu ");
        next_cell_info = next_cell_info2;

        mdl.addCell(next_cell_info.1,sp_heat_cap,density,kx,ky,kz,init_temp,
                    next_cell_info.3,next_cell_info.2,next_cell_info.4, &mut temporary_templist,
                    next_cell_info.5,next_cell_info.6, &mut cell_indices);
        // global_time = global_time + time;
        let timedt = 0.1;
        for i in 0..(60.0/timedt) as usize{
            global_time = global_time + timedt;
            mdl.run_model(timedt,dt,conv_coeff,t_env,bw, &mut global_time,&pool,areas_and_dists, &mut temporary_templist,
                          conductivity_interp,maxthreads,&cell_indices);

        }
        println!("{} out of {}, time {}",endval, endval, global_time);
    }

}
