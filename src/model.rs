use crate::primitives::Node;
use std::collections::HashMap;
use crate::interpolator::Interpolator;
use rayon::ThreadPool;
use std::io::BufWriter;
use std::fs::File;
use crate::model_updater::ModelUpdater;
use std::io::Write;
use std::cmp::max;

/// Model has list of Nodes, list of cells, list of temperature for each cell, and information about
/// six neighbors (other cell or free surface) for each cell
pub struct Model{
    nodelist:Vec<Node>,
    sp_heat_cap:Vec<f64>,
    density:Vec<f64>,
    kx:Vec<f64>,
    ky:Vec<f64>,
    kxy:Vec<f64>,
    kz:Vec<f64>,
    pub(crate) templist:Vec<f64>,
    pub(crate) templist_old:Vec<f64>,
    has_neighbor:Vec<[bool;6]>,
    neighbors:Vec<[isize;6]>,
    pub(crate) temp_at_contact:Vec<[f64;6]>,
    density_vol_sp_ht_cap:Vec<f64>,
    pub(crate) cell_index_map:HashMap<isize,isize>,
    sp_heat_cap_interp:Interpolator,
    beadwidth:f64,
    beadheight:f64,
    active_layer_first_element:usize,
    previous_active_layer:usize,
    previous_layer:usize,
    current_layer:usize,
    turn_off_layers_at:usize,
    first_element_of_each_layer:Vec<usize>,
    temp_pos_for_split:Vec<usize>,
    pub(crate)interface: Vec<[usize;2]>,
    pub(crate)interface_interlayer:Vec<bool>,
    pub(crate)interface_interbead:Vec<bool>,
    pub(crate)interface_within_bead:Vec<bool>,
    pub (crate) degree_bonding: Vec<f64>,
    tr_interpolation_table: Interpolator,
    interface_temps:Vec<f64>,
}

/// methods for Model struct
impl Model{
    /// generates a new empty Model
    pub fn new(nodes:Vec<Node>,activation_cells_length:usize, sp_heat_cap_interp: Interpolator, beadwidth:f64, beadheight:f64,
               turn_off_layers_at:usize, tr_interpolation_table:Interpolator)->Model{
        return Model{
            nodelist: nodes,
            sp_heat_cap:Vec::with_capacity(activation_cells_length),
            density:Vec::with_capacity(activation_cells_length),
            kx:Vec::with_capacity(activation_cells_length),
            ky:Vec::with_capacity(activation_cells_length),
            kxy:Vec::with_capacity(activation_cells_length),
            kz:Vec::with_capacity(activation_cells_length),
            templist: Vec::with_capacity(activation_cells_length),
            templist_old:Vec::with_capacity(activation_cells_length),
            neighbors: Vec::with_capacity(activation_cells_length),
            has_neighbor:Vec::with_capacity(activation_cells_length),
            temp_at_contact: Vec::with_capacity(activation_cells_length),
            density_vol_sp_ht_cap:Vec::with_capacity(activation_cells_length),
            //celllist:Vec::with_capacity(activation_cells_length);
            cell_index_map:HashMap::with_capacity(activation_cells_length),
            sp_heat_cap_interp,
            beadwidth,
            beadheight,
            active_layer_first_element:0,
            previous_active_layer:0,
            previous_layer: 0,
            current_layer: 0,
            turn_off_layers_at,
            first_element_of_each_layer: vec![0],
            temp_pos_for_split:Vec::with_capacity(activation_cells_length),
            interface: Vec::with_capacity(activation_cells_length*4),
            degree_bonding: Vec::with_capacity(activation_cells_length*4),
            tr_interpolation_table: tr_interpolation_table,
            interface_temps:Vec::with_capacity(activation_cells_length*4),
            interface_interlayer:Vec::with_capacity(activation_cells_length*4),
            interface_interbead:Vec::with_capacity(activation_cells_length*4),
            interface_within_bead:Vec::with_capacity(activation_cells_length*4),
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
        if (layer - self.previous_active_layer)>self.turn_off_layers_at{
            self.previous_active_layer = self.previous_active_layer+1;
            self.active_layer_first_element = self.first_element_of_each_layer[self.previous_active_layer+1];

        }

        let p1 = self.nodelist[node_no[0]].clone();
        let p2 = self.nodelist[node_no[7]].clone();

        let volume = ((p2.pt.x-p1.pt.x)*(p2.pt.y-p1.pt.y)*(p2.pt.z-p1.pt.z)).abs();

        self.density_vol_sp_ht_cap.push(density*volume);
        temporary_templist.push(init_temp.clone());


        // please check this one more time
        self.kx.push(kx);
        self.ky.push(ky);//+kx*costheta*sintheta-ky*costheta*sintheta);
        //self.kxy.push((kx-ky)*costheta*sintheta);
        //
        self.kz.push(kz);
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
        /* let mut is_neighbor_f64:[f64;6] = [0.0,0.0,0.0,0.0,0.0,0.0];
         let mut is_no_neighbor_f64 = [1.0,1.0,0.1,0.0,1.0,1.0];
         for i in 0..6{
             if is_neighbor[i]{
                 is_neighbor_f64[i] = 1.0;
                 is_no_neighbor_f64[i] = 0.0;
             }
             else{
                 is_neighbor_f64[i] = 0.0;
                 is_no_neighbor_f64[i] = 1.0;
             }
         }*/
        // println!("{:?}",neighbor2);
        // println!("{:?}",is_neighbor);
        self.templist.push(init_temp);
        self.templist_old.push(init_temp);
        if is_neighbor[0] {
            self.neighbors[neighbor2[0] as usize][1] = counter as isize;
            self.has_neighbor[neighbor2[0] as usize][1] = true;
            self.temp_at_contact[neighbor2[0] as usize][1] = self.templist[neighbor2[0] as usize];
            self.interface.push([counter,neighbor2[0] as usize]);
            self.interface_temps.push(0.0_f64);
            self.degree_bonding.push(0.0_f64);
            self.interface_within_bead.push(true);
            self.interface_interlayer.push(false);
            self.interface_interbead.push(false);
            //self.temp_at_contact[neighbor2[0] as usize][1] = self.templist[neighbor2[0] as usize][1]
            //  self.has_neighbor_f64[neighbor2[0] as usize][1] = 1.0;
            //  self.has_no_neighbor_f64[neighbor2[0] as usize][1] = 0.0;
        }
        //println!("check 1 pass");
        if is_neighbor[1] {
            self.neighbors[neighbor2[1] as usize][0]=counter as isize;
            self.temp_at_contact[neighbor2[1] as usize][0] = self.templist[neighbor2[1] as usize];
            self.has_neighbor[neighbor2[1] as usize][0]=true;
            self.interface.push([counter,neighbor2[1] as usize]);
            self.interface_temps.push(0.0_f64);
            self.degree_bonding.push(0.0_f64);
            self.interface_within_bead.push(true);
            self.interface_interlayer.push(false);
            self.interface_interbead.push(false);
            // self.has_neighbor_f64[neighbor2[1] as usize][0] = 1.0;
            //self.has_no_neighbor_f64[neighbor2[1] as usize][0] = 0.0;
        }
        //  println!("check 2 pass");
        if is_neighbor[2] {
            self.neighbors[neighbor2[2] as usize][3]=counter as isize;
            self.temp_at_contact[neighbor2[2] as usize][3] = self.templist[neighbor2[2] as usize];
            self.has_neighbor[neighbor2[2] as usize][3]=true;
            self.interface.push([counter,neighbor2[2] as usize]);
            self.interface_temps.push(0.0_f64);
            self.degree_bonding.push(0.0_f64);
            self.interface_within_bead.push(false);
            self.interface_interlayer.push(false);
            self.interface_interbead.push(true);
            // self.has_neighbor_f64[neighbor2[2] as usize][3] = 1.0;
            // self.has_neighbor_f64[neighbor2[2] as usize][3] = 0.0;
        }
        // println!("check 3 pass");
        if is_neighbor[3] {
            self.neighbors[neighbor2[3] as usize][2]=counter as isize;
            self.temp_at_contact[neighbor2[3] as usize][2] = self.templist[neighbor2[3] as usize];
            self.has_neighbor[neighbor2[3] as usize][2] = true;
            self.interface.push([counter,neighbor2[3] as usize]);
            self.interface_temps.push(0.0_f64);
            self.degree_bonding.push(0.0_f64);
            self.interface_within_bead.push(false);
            self.interface_interlayer.push(false);
            self.interface_interbead.push(true);
            //self.has_neighbor_f64[neighbor2[3] as usize][2] = 1.0;
            //self.has_neighbor_f64[neighbor2[3] as usize][2] = 0.0;
        }
        // println!("check 4 pass");
        if is_neighbor[4] {
            self.neighbors[neighbor2[4] as usize][5]=counter as isize;
            self.temp_at_contact[neighbor2[4] as usize][5] = self.templist[neighbor2[4] as usize];
            self.has_neighbor[neighbor2[4] as usize][5]=true;
            self.interface.push([counter,neighbor2[4] as usize]);
            self.interface_temps.push(0.0_f64);
            self.degree_bonding.push(0.0_f64);
            self.interface_within_bead.push(false);
            self.interface_interlayer.push(true);
            self.interface_interbead.push(false);
        }

        if is_neighbor[5] {
            self.neighbors[neighbor2[5] as usize][4]= counter as isize;
            self.temp_at_contact[neighbor2[5] as usize][4] = self.templist[neighbor2[5] as usize];
            self.has_neighbor[neighbor2[5] as usize][4] = true;
            self.interface.push([counter,neighbor2[5] as usize]);
            self.interface_temps.push(0.0_f64);
            self.degree_bonding.push(0.0_f64);
            self.interface_within_bead.push(false);
            self.interface_interlayer.push(true);
            self.interface_interbead.push(false);
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
        }
    }




    /// calculates the temperature of the cells in the model during a time period with a given time
    /// step. Takes time period, time step, and a file (buffered) where output is to be written as
    /// input. Writes output to the file (buffer). Might need to change that for speed.
    pub fn run_model(&mut self, time:f64, dt:f64, conv_coeff:f64, t_env:f64, file:&mut BufWriter<File>,
                     global_time:&mut f64,pool:&ThreadPool, areas_and_dists:[f64;5], newtemp: &mut Vec<f64>,
                     maxthreads:usize, cell_indices:&Vec<usize>, datastorer:&mut (Vec<Vec<f64>>,Vec<Vec<f64>>) ) {

        let count = (time / dt) as usize;
        let (inactive_temps,active_temps) = newtemp.split_at_mut(self.active_layer_first_element);
        let (old_cell_indices_slice, new_cell_indices_slice) = cell_indices.split_at(self.active_layer_first_element);
        let new_cell_slice_start = new_cell_indices_slice[0].clone();
        let ncis_len = new_cell_indices_slice.len();
        let new_cell_slice_end = new_cell_indices_slice[ncis_len-1].clone();

        for i in 0..count {
            self.find_new_cell_temp_all(dt, conv_coeff, t_env, areas_and_dists[0], areas_and_dists[1], areas_and_dists[2],
                                        areas_and_dists[3], areas_and_dists[4], &pool,
                                        active_temps,1,maxthreads, new_cell_indices_slice);


            self.templist[new_cell_slice_start..new_cell_slice_end+1].copy_from_slice(active_temps);
            // Model::parallel_copy_newtemps(active_temps,&mut self.templist[new_cell_slice_start..new_cell_slice_end+1],1,maxthreads,pool);
            //if (self.templist.len()>0) {write!(file, "{},{},{},{}\n", self.templist[0].clone(), global_time,self.sp_heat_cap[0],self.density_vol_sp_ht_cap[0] );}


        }
        self.find_new_cell_temp_all(time - count as f64 * dt, conv_coeff, t_env, areas_and_dists[0], areas_and_dists[1], areas_and_dists[2],
                                    areas_and_dists[3], areas_and_dists[4], &pool,
                                    active_temps,1,maxthreads, new_cell_indices_slice);
        //self.templist_old[new_cell_slice_start..new_cell_slice_end+1].copy_from_slice(&self.templist[new_cell_slice_start..new_cell_slice_end+1]);

        // let zz = 0.0 as f64;
        /*for i in new_cell_slice_start..new_cell_slice_end+1{
           //let kk = active_temps[i-new_cell_slice_start];
            //let kk2 = self.templist_old[i];
           // let diffkk = (kk-kk2).abs();
           // let check = diffkk > 0.5;
            if (active_temps[i-new_cell_slice_start]-self.templist_old[i]).abs()>0.5{
                self.templist_old[i]  = active_temps[i-new_cell_slice_start];
                (*datastorer).0[i].push(*global_time);
                (*datastorer).1[i].push(self.templist[i]);
            }
        }*/

        //Model::parallel_write_datastorer(active_temps, &mut self.templist_old[new_cell_slice_start..new_cell_slice_end+1],
        // &mut self.templist[new_cell_slice_start..new_cell_slice_end+1], 1, maxthreads,
        //                               datastorer.0.as_mut_slice(),datastorer.1.as_mut_slice(), pool, global_time);




        // update degree of healing
        // find interface temp
        //println!("{}",self.sp_heat_cap.len());
        for i in 0..self.interface_temps.len(){
            let c1 = self.interface[i][0];
            let c2 = self.interface[i][1];
            let interface_temp = (self.templist[c1]+self.templist[c2])/2.0_f64;
            if (interface_temp > self.tr_interpolation_table.xvalues[0]) {
                //let kk = (interface_temp-self.tr_interpolation_table.xvalues[0])/self.tr_interpolation_table.xstep;
                //let kk2 = interface_temp-self.tr_interpolation_table.xvalues[0];
                let location = ((interface_temp - self.tr_interpolation_table.xvalues[0]) / self.tr_interpolation_table.xstep) as usize;
                let x1 = self.tr_interpolation_table.xvalues[location];
                let x2 = self.tr_interpolation_table.xvalues[location + 1];
                let y1 = self.tr_interpolation_table.yvalues[location];
                let y2 = self.tr_interpolation_table.yvalues[location + 1];
                let slope = (y2 - y1) / (x2 - x1);
                let intercept = y1 - slope * x1;
                let rep_time = slope * interface_temp + intercept;
                let db = (time / rep_time).sqrt();
                self.degree_bonding[i] = self.degree_bonding[i] + db;
            }

        }

        for i in 0..self.degree_bonding.len(){
            if self.degree_bonding[i] > 1.0_f64{
                self.degree_bonding[i] = 1.0_f64;
            }
        }
        // find reptation time for each interface at its current temp
        self.templist[new_cell_slice_start..new_cell_slice_end+1].copy_from_slice(active_temps);




        // (*datastorer).0.push(global_time.clone());
        //let templist_clone = self.templist.clone().as_slice();
        //(*datastorer).1.push(active_temps.to_vec());
        //(*datastorer).2.push(new_cell_slice_start.clone());
    }

    /// parallel copy to newtemps from active temps
    fn parallel_copy_newtemps(active_temps:&[f64], newtemps:&mut[f64], num_threads:usize, max_threads:usize, pool:&ThreadPool){
        if max_threads < num_threads{
            let templen = active_temps.len();
            let (at0, at1) = active_temps.split_at(templen/2);
            let (nt0, nt1) = newtemps.split_at_mut(templen/2);
            pool.install(||rayon::join(
                || Model::parallel_copy_newtemps(at0,nt0,num_threads*2, max_threads, pool),
                || (Model::parallel_copy_newtemps(at1,nt1,num_threads*2, max_threads,pool))));
        }
        else{
            newtemps.copy_from_slice(active_temps);
        }

    }
    ///parallel write to datastorer
    fn parallel_write_datastorer(active_temps:&[f64], oldtemps:&mut[f64],newtemps:&mut[f64], num_threads:usize, max_threads:usize,
                                 ds0:&mut [Vec<f64>], ds1:&mut [Vec<f64>], pool: &ThreadPool, globaltime:&f64)
    {
        if max_threads < num_threads{
            let dslen = ds0.len();
            let (ds00, ds01) = ds0.split_at_mut(dslen/2);
            let (ds10, ds11) = ds1.split_at_mut(dslen/2);
            let (ot0,ot1) = oldtemps.split_at_mut(dslen/2);
            let (nt0,nt1) = newtemps.split_at_mut(dslen/2);
            let(at0, at1) = active_temps.split_at(dslen/2);
            pool.install(||rayon::join(
                || Model::parallel_write_datastorer( at0, ot0, nt0, num_threads*2, max_threads, ds00, ds10, pool, globaltime),
                || Model::parallel_write_datastorer(at1, ot1, nt1, num_threads*2, max_threads, ds01, ds11, pool, globaltime)
            ));
        }
        else{
            for i in 0..active_temps.len(){
                if (active_temps[i]-oldtemps[i]).abs() > 0.5 {
                    oldtemps[i] = active_temps[i];
                    ds0[i].push(*globaltime);
                    ds1[i].push(active_temps[i]);
                }
            }
            newtemps.copy_from_slice(active_temps);
        }

    }

    /// updates model based on input. Needs to be placed in Model implementation later.
    pub fn update_model(mu:&mut ModelUpdater, mdl:&mut Model, bw:&mut BufWriter<File>, areas_and_dists:[f64;5], time_step:f64,
                        maxthreads:usize, input_data:([f64;3],f64,f64,f64,[f64;2]), pool: &ThreadPool) -> (Vec<Vec<f64>>, Vec<Vec<f64>>)
    {


        //let pool = rayon::ThreadPoolBuilder::new().num_threads(maxthreads).build().unwrap();
        let mut temporary_templist = Vec::with_capacity(mu.activation_times.len());
        let mut cell_indices = Vec::with_capacity(mu.activation_times.len());
        let mut store_temp_data = Vec::with_capacity(mu.activation_times.len());
        let mut store_global_time = Vec::with_capacity(mu.activation_times.len());
        //let mut store_starting_cell = Vec::with_capacity(mu.activation_times.len());
        let mut datastorer = (store_global_time, store_temp_data);
        println!("activation times len {}" ,mu.activation_times.len());
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
                    next_cell_info.3,next_cell_info.2,next_cell_info.4, &mut temporary_templist,next_cell_info.5,next_cell_info.6, &mut cell_indices);
        // let mut cts_dat = Vec::with_capacity(100000);
        // cts_dat.push([global_time,init_temp]);
        // mdl.temp_data_store.push(cts_dat);
        // global_time = global_time+time;
        datastorer.0.push(Vec::with_capacity(1000));
        datastorer.1.push(Vec::with_capacity(1000));
        mdl.run_model(time,dt,conv_coeff,t_env,bw, &mut global_time,&pool,areas_and_dists, &mut temporary_templist,
                      maxthreads,&cell_indices, &mut datastorer);

        //println!("run model");
        let endval = mu.activation_times.len();
        for i in 1..endval-1 {

            next_cell_info = next_cell_info2;
            next_cell_info2 = mu.get_next_cell_info();

            let time = next_cell_info2.0 - next_cell_info.0;
            // global_time = global_time + time;
            if i%1000 == 0 {println!("{} out of {}, timeperiod {}",i, endval, time);}
            //println!("layer {}",next_cell_info.6);
            mdl.addCell(next_cell_info.1, sp_heat_cap, density, kx, ky, kz,init_temp,
                        next_cell_info.3, next_cell_info.2,next_cell_info.4, &mut temporary_templist,
                        next_cell_info.5,next_cell_info.6, &mut cell_indices);
            datastorer.0.push(Vec::with_capacity(2000));
            datastorer.1.push(Vec::with_capacity(2000));
            //let mut cts_dat = Vec::with_capacity(100000);
            //cts_dat.push([global_time,init_temp]);
            //mdl.temp_data_store.push(cts_dat);
            //println!("second addcell");
            //println!("second run_model");
            global_time = global_time+time;
            mdl.run_model(time, dt, conv_coeff, t_env,  bw, &mut global_time,&pool,areas_and_dists,
                          &mut temporary_templist, maxthreads,&cell_indices, &mut datastorer);

            //println!("second run model");
        }
        //println!("pulupulu ");
        next_cell_info = next_cell_info2;

        mdl.addCell(next_cell_info.1, sp_heat_cap, density, kx, ky, kz,init_temp,
                    next_cell_info.3, next_cell_info.2,next_cell_info.4, &mut temporary_templist,
                    next_cell_info.5,next_cell_info.6,&mut cell_indices);
        datastorer.0.push(Vec::with_capacity(2000));
        datastorer.1.push(Vec::with_capacity(2000));
        ///let mut cts_dat = Vec::with_capacity(100000);
        //cts_dat.push([global_time,init_temp]);
        //mdl.temp_data_store.push(cts_dat);
        // global_time = global_time + time;
        for i in 0..60 {
            global_time = global_time+1.0;
            mdl.run_model(1.0, dt, conv_coeff, t_env, bw, &mut global_time, &pool, areas_and_dists, &mut temporary_templist,
                          maxthreads,&cell_indices, &mut datastorer);

        }
        println!("{} out of {}, timeperiod {}",endval, endval, time);
        return datastorer;
    }

}
