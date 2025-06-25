use crate::primitives::Node;
use std::collections::HashMap;
use crate::interpolator::Interpolator;
use rayon::ThreadPool;
use std::io::BufWriter;
use std::fs::File;
use crate::model_updater::ModelUpdater;
use std::io::Write;
use std::borrow::{BorrowMut, Borrow};


const STEFAN_BOLTZMAN:f64 = 5.67e-8;

/// Model has list of Nodes, list of cells, list of temperature for each cell, and information about
/// six neighbors (other cell or free surface) for each cell
pub struct Model{
    nodelist:Vec<[f64;3]>,
    nodetemp: Vec<f64>,
    nodetemp_old: Vec<f64>,
    total_nodes_activated:usize,
    node_activated: Vec<bool>,
    is_bot_cell: Vec<bool>,
    bed_k: f64,
    bed_temp:f64,
    sp_heat_cap:Vec<f64>,
    density:Vec<f64>,
    kx:Vec<f64>,
    ky:Vec<f64>,
    kz:Vec<f64>,
    h: Vec<f64>,
    pub(crate) templist:Vec<f64>,
    pub(crate) oldtemplist:Vec<f64>,
    has_neighbor:Vec<[bool;6]>,
    neighbors:Vec<[usize;6]>,
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
    element_no:Vec<usize>,
    zmin:f64,
    emissivity:f64,

}

/// methods for Model struct
impl Model{
    /// generates a new empty Model
    pub fn new(nodes:Vec<[f64;3]>,activation_cells_length:usize, sp_heat_cap_interp: Interpolator, element_width:f64, element_height:f64,
               turn_off_layers_at:usize, zmin: f64, density:f64, emissivity:f64, bed_k:f64, bed_temp:f64, maxthreads: usize)->Model{
        println!("element width {} element height {}", element_width, element_height);
        let nodeslen = nodes.len();
        return Model{
            nodelist: nodes,
            nodetemp: Vec::with_capacity(nodeslen),
            nodetemp_old: Vec::with_capacity(nodeslen),
            node_activated: vec![false; nodeslen],
            is_bot_cell: vec![false;activation_cells_length],
            bed_k,
            bed_temp,
            total_nodes_activated:0,
            sp_heat_cap:Vec::with_capacity(activation_cells_length),
            density:Vec::with_capacity(activation_cells_length),
            kx:Vec::with_capacity(activation_cells_length),
            ky:Vec::with_capacity(activation_cells_length),
            kz:Vec::with_capacity(activation_cells_length),
            h:Vec::with_capacity(activation_cells_length),
            templist: Vec::with_capacity(activation_cells_length),
            oldtemplist: Vec::with_capacity(activation_cells_length),
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
            temp_pos_for_split: Vec::with_capacity(activation_cells_length),
            emissivity:emissivity,
            // strain_time_dep: Vec::with_capacity(activation_cells_length),
            // strain_time_dep_prev: Vec::with_capacity(activation_cells_length),
            // stress: Vec::with_capacity(activation_cells_length),
            // stress_prev: Vec::with_capacity(activation_cells_length),
            // loads: Vec::with_capacity(activation_cells_length),
            element_no:Vec::with_capacity(activation_cells_length),
            // stm: StructuralModel::with_capacity(activation_cells_length,element_width,element_width,element_height, density),
            // gbcs: Vec::with_capacity(activation_cells_length),
            zmin,
            //  Evec: Vec::with_capacity(activation_cells_length),
            // tmp1,
            // tmp2,
            // Ctmp,
        }
    }



    /// adds a cell to the model. Takes node index, specific heat capacity, density, conductivity in
    /// the x-direction, conductivity in the y and z-directions, initial temperature, and neighbor
    /// information of the cell
    pub fn addCell(&mut self,node_no:[usize;8],sp_heat_cap:f64,
                   density:f64, kx:f64, ky:f64, kz:f64, init_temp:f64,is_neighbor:[bool;6],
                   neighbor:[isize;6],cell_index:usize,temporary_templist:&mut Vec<f64>,
                   orientation:[f64;2],layer:usize,cell_indices:&mut Vec<usize>, node_data_storer: &mut Vec<Vec<f64>>)
    {
        let mut counter = self.sp_heat_cap.len();
        self.element_no.push(counter);
        self.cell_index_map.insert(cell_index as isize,counter as isize);
        self.templist.push(init_temp);
        
        let mut cellz = 0.0;
        for i in 0..8{
            cellz += self.nodelist[node_no[i]][2];
        }
        cellz = cellz / 8.0;
        //variable h
        // h = 3 at z = 0;
        // h = 15 at z = 1;
        
        self.h.push(3.0 + (cellz - 0.0) / (1.0 - 0.0) * (cellz - 0.0) * (15.0 - 3.0));

        self.oldtemplist.push(0.0);
        temporary_templist.push(init_temp);
        if layer == 0{
            //println!("layer 0 element");
            self.is_bot_cell[counter] = true;
        }

        //self.Evec.push(1e8);
        // let mut bcs = [false;24];
        for kk in 0..8 {

            let zval = self.nodelist[node_no[kk]][2];
            if !self.node_activated[node_no[kk]]{
                self.node_activated[node_no[kk]] = true;
                self.total_nodes_activated += 1;
                node_data_storer.push(Vec::with_capacity(100));
                self.nodetemp.push(init_temp);
                self.nodetemp_old.push(0.0);
            }
            //if (zval-self.zmin).abs() < 1e-7{
            //   bcs[kk*3] = true;
            //    bcs[kk*3+1] = true;
            //    bcs[kk*3+2] = true;
            // }
        }
        let elem = [node_no[0],node_no[1],node_no[3],node_no[2],node_no[4],node_no[5],node_no[7],node_no[6]];
        /*println!("zmin, {}", self.zmin);
        println!("element_coord start");
        if self.templist.len()<27{
            for i in 0..8{
                println!("{},{},{}", &self.nodelist[elem[i]].pt.x,&self.nodelist[elem[i]].pt.y,&self.nodelist[elem[i]].pt.z);
                println!("{},{},{}",&bcs[i*3],&bcs[i*3+1],&bcs[i*3+2]);
            }
        }
        println!("element_coord end");*/



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
        // temporary_templist.push(init_temp.clone());
        // temporary_Glist.push(0.0);
        //temporary_nlist.push(0.0);
        // temporary_tdlist.push(0.0);
        // self.strain_time_dep.push(0.0);
        // self.strain_time_dep_prev.push(0.0);
        // self.stress.push(0.0);
        //  self.stress_prev.push(0.0);



        // please check this one more time
        self.kx.push(kx);
        self.ky.push(ky);
        self.kz.push(kz);

        //self.kz.push(kz);
        self.density.push(density);
        self.sp_heat_cap.push(sp_heat_cap);
        // println!("start suspicious code");

        //update neighbor index
        let mut neighbor2:[usize;6] =[0,0,0,0,0,0];
        self.temp_at_contact.push([-100.0;6]);
        for i in 0..6{
            if is_neighbor[i]{
                neighbor2[i] = self.cell_index_map.get(&(neighbor[i])).expect("can't find it").clone() as usize;
            }
        }

        //self.templist.push(init_temp);
        if is_neighbor[0] {
            self.neighbors[neighbor2[0] as usize][1] = counter;
            self.has_neighbor[neighbor2[0] as usize][1] = true;
            self.temp_at_contact[neighbor2[0] as usize][1] = self.templist[neighbor2[0] as usize];
            //self.temp_at_contact[neighbor2[0] as usize][1] = self.templist[neighbor2[0] as usize][1]
            //  self.has_neighbor_f64[neighbor2[0] as usize][1] = 1.0;
            //  self.has_no_neighbor_f64[neighbor2[0] as usize][1] = 0.0;
        }
        //println!("check 1 pass");
        if is_neighbor[1] {
            self.neighbors[neighbor2[1] as usize][0]=counter;
            self.temp_at_contact[neighbor2[1] as usize][0] = self.templist[neighbor2[1] as usize];
            self.has_neighbor[neighbor2[1] as usize][0]=true;
            // self.has_neighbor_f64[neighbor2[1] as usize][0] = 1.0;
            //self.has_no_neighbor_f64[neighbor2[1] as usize][0] = 0.0;
        }
        //  println!("check 2 pass");
        if is_neighbor[2] {
            self.neighbors[neighbor2[2] as usize][3]=counter ;
            self.temp_at_contact[neighbor2[2] as usize][3] = self.templist[neighbor2[2] as usize];
            self.has_neighbor[neighbor2[2] as usize][3]=true;
            // self.has_neighbor_f64[neighbor2[2] as usize][3] = 1.0;
            // self.has_neighbor_f64[neighbor2[2] as usize][3] = 0.0;
        }
        // println!("check 3 pass");
        if is_neighbor[3] {
            self.neighbors[neighbor2[3] as usize][2]=counter ;
            self.temp_at_contact[neighbor2[3] as usize][2] = self.templist[neighbor2[3] as usize];
            self.has_neighbor[neighbor2[3] as usize][2] = true;
            //self.has_neighbor_f64[neighbor2[3] as usize][2] = 1.0;
            //self.has_neighbor_f64[neighbor2[3] as usize][2] = 0.0;
        }
        // println!("check 4 pass");
        if is_neighbor[4] {
            self.neighbors[neighbor2[4] as usize][5]=counter;
            self.temp_at_contact[neighbor2[4] as usize][5] = self.templist[neighbor2[4] as usize];
            self.has_neighbor[neighbor2[4] as usize][5]=true;
        }

        if is_neighbor[5] {
            self.neighbors[neighbor2[5] as usize][4]= counter;
            self.temp_at_contact[neighbor2[5] as usize][4] = self.templist[neighbor2[5] as usize];
            self.has_neighbor[neighbor2[5] as usize][4] = true;
        }

        self.has_neighbor.push(is_neighbor);
        self.neighbors.push(neighbor2);
        //self.stm.add_element(&self.nodelist,elem, &bcs, self.has_neighbor[counter][5]);


        //go through all the active layers, add weights to cells directly below this cell being deposited
        // self.loads.push(-volume * &self.density[counter] * 9.81); //9.81 is g
        //self.loads.push(0.0);
        //let mut endcell = false;
        //let mut cellcount = 0;
        // let mut botcell:usize = 0;
        /* if !self.has_neighbor[counter][5] {
             endcell=true
         }else{
             botcell = self.neighbors[counter][5] as usize;
         }
         while !endcell && cellcount < self.turn_off_layers_at{
             self.loads[botcell] -= volume * self.density[botcell] * 9.81; //g is 9.81

             if !self.has_neighbor[botcell][5] {
                 endcell=true
             }else{
                 botcell = self.neighbors[botcell][5] as usize;
             }
         }
         /*println!("layer {}",self.current_layer);
         if self.loads.len() <1000 && self.current_layer != self.previous_layer {
             for i in 0..self.loads.len() {
                 println!("load,{},{}", i, self.loads[i]);
             }
         }*/*/

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
    pub fn find_new_temp_all_split(&self, dt:f64, conv_coeff2:f64, t_env:f64, temp: &mut [f64], side_area:f64, top_area:f64, sa_d_sd:f64, ta_d_td:f64, vol:f64,
                                   cell_indices_slice:&[usize]){
        for cell in cell_indices_slice{
            let mut dQ = 0.0;
            let thiscelltemp = self.templist[*cell];
            
           
            // T^4 needed for radiation calculations
            let T4 = thiscelltemp*thiscelltemp*thiscelltemp*thiscelltemp;
            let Ta4 = t_env*t_env*t_env*t_env;
            //println!("dt conv_coeff t_env thiscelltemp top_area side_area sa_d_sd ta_d_td vol {} {} {} {} {} {} {} {} ", dt, conv_coeff, thiscelltemp, top_area, side_area, sa_d_sd, ta_d_td, self.density_vol_sp_ht_cap[*cell]);
            //front back
            //dist area and vol is  [ side area, top area, sidearea/sidedistance, toparea/topdistance, volume]
            unsafe {
                let conv_coeff = self.h.get_unchecked(*cell);
                for i in 0..2 {
                    if !self.has_neighbor.get_unchecked(*cell).get_unchecked(i) {
                        dQ = dQ - (conv_coeff * (self.templist.get_unchecked(*cell) - t_env) + self.emissivity * (T4 - Ta4) * STEFAN_BOLTZMAN) * dt * side_area;// * self.has_no_neighbor_f64[*cell][i] ;
                    } else {
                        let neighbor_index = self.neighbors.get_unchecked(*cell).get_unchecked(i);
                        let t1 = self.templist.get_unchecked(*cell);
                        let t2 = self.templist.get_unchecked(*neighbor_index);
                        dQ = dQ - self.kx.get_unchecked(*cell) * (t1 - t2) * sa_d_sd * dt;// *self.has_neighbor_f64[*cell][i];
                    }
                }
                //right left
                for i in 2..4 {
                    if !self.has_neighbor.get_unchecked(*cell).get_unchecked(i) {
                        dQ = dQ - (conv_coeff * (self.templist.get_unchecked(*cell) - t_env) + self.emissivity * (T4 - Ta4) * STEFAN_BOLTZMAN) * dt * side_area;// * self.has_no_neighbor_f64[*cell][i] ;
                    } else {
                        let neighbor_index = self.neighbors.get_unchecked(*cell).get_unchecked(i);
                        let t1 = self.templist.get_unchecked(*cell);
                        let t2 = self.templist.get_unchecked(*neighbor_index);
                        dQ = dQ - self.ky.get_unchecked(*cell) * (t1 - t2) * sa_d_sd * dt;// *self.has_neighbor_f64[*cell][i];
                    }
                }

                //top down
                for i in 4..5 {
                    if !self.has_neighbor.get_unchecked(*cell).get_unchecked(i) {
                        dQ = dQ - (conv_coeff * (self.templist.get_unchecked(*cell) - t_env) + self.emissivity * (T4 - Ta4) * STEFAN_BOLTZMAN) * dt * top_area;
                    } else {
                        let neighbor_index = self.neighbors.get_unchecked(*cell).get_unchecked(i);
                        let t1 = self.templist.get_unchecked(*cell);
                        let t2 = self.templist.get_unchecked(*neighbor_index);
                        dQ = dQ - self.kz.get_unchecked(*cell) * (t1 - t2) * ta_d_td * dt;
                    }
                }
                // bottom layer code check conduction to bed
                for i in 5..6 {
                    if *self.is_bot_cell.get_unchecked(*cell) {
                        dQ = dQ - self.bed_k * (self.templist.get_unchecked(*cell) - self.bed_temp) * ta_d_td * dt;
                    } else if !self.has_neighbor.get_unchecked(*cell).get_unchecked(i) {
                        dQ = dQ - (conv_coeff * (self.templist.get_unchecked(*cell) - t_env) + self.emissivity * (T4 - Ta4) * STEFAN_BOLTZMAN) * dt * top_area;
                    } else {
                        let neighbor_index = self.neighbors.get_unchecked(*cell).get_unchecked(i);
                        let t1 = self.templist.get_unchecked(*cell);
                        let t2 = self.templist.get_unchecked(*neighbor_index);
                        dQ = dQ - self.kz.get_unchecked(*cell) * (t1 - t2) * ta_d_td * dt;
                    }
                }

                let mult = self.density_vol_sp_ht_cap.get_unchecked(*cell) * self.sp_heat_cap.get_unchecked(*cell);
                //println!("cell indices slice {:?}", cell_indices_slice);

                *temp.get_unchecked_mut(*cell - cell_indices_slice.get_unchecked(0)) = (mult * self.templist.get_unchecked(*cell) + dQ) / mult;
            }

        }

    }


    /// calculates the temperature of the cells in the model during a time period with a given time
    /// step. Takes time period, time step, and a file (buffered) where output is to be written as
    /// input. Writes output to the file (buffer). Might need to change that for speed.
    pub fn run_model(&mut self, time:f64, dt1:f64, conv_coeff:f64, t_env:f64, file:&mut BufWriter<File>,
                     global_time:&mut f64,pool:&ThreadPool, areas_and_dists:[f64;5], newtemp: &mut Vec<f64>,
                     maxthreads:usize, cell_indices:&Vec<usize>
    ) {
        let mut dt = dt1;
        if time<dt1 {dt = time;}

        let count = (time / dt) as usize;


        let (inactive_temps, active_temps) = newtemp.split_at_mut(self.active_layer_first_element);
        let (old_cell_indices_slice, new_cell_indices_slice) = cell_indices.split_at(self.active_layer_first_element);
        let new_cell_slice_start = new_cell_indices_slice[0].clone();
        let ncis_len = new_cell_indices_slice.len();
        let new_cell_slice_end = new_cell_indices_slice[ncis_len - 1].clone();
        {
            let (inactive_selftemps, active_selftemps) = self.templist.split_at(self.active_layer_first_element);
            let (sp_ht_cap_inactive, sp_ht_cap_active) = self.sp_heat_cap.split_at_mut(self.active_layer_first_element);

            self.sp_heat_cap_interp.interpolate_inner(active_selftemps, sp_ht_cap_active, 1, maxthreads, pool);

            let (inactive_conductivity, active_conductivity) = self.kx.split_at_mut(self.active_layer_first_element);
        }

        //println!("time {}",time);
        for i in 0..count {
            self.find_new_cell_temp_all(dt, conv_coeff, t_env, areas_and_dists[0], areas_and_dists[1], areas_and_dists[2],
                                        areas_and_dists[3], areas_and_dists[4], &pool,
                                        active_temps, 1, maxthreads, new_cell_indices_slice);

            self.templist[new_cell_slice_start..new_cell_slice_end + 1].copy_from_slice(active_temps);
        }
        self.find_new_cell_temp_all(time - count as f64 * dt, conv_coeff, t_env, areas_and_dists[0], areas_and_dists[1], areas_and_dists[2],
                                    areas_and_dists[3], areas_and_dists[4], &pool,
                                    active_temps, 1, maxthreads, new_cell_indices_slice);

    }




    fn parallel_copy(valold: &mut [f64], valnew: &[f64], numthreads: usize, maxthreads: usize, pool: &ThreadPool)
    {
        if numthreads < maxthreads {
            let splitpos = valold.len() / 2;
            let (o1, o2) = valold.split_at_mut(splitpos);
            let (n1, n2) = valnew.split_at(splitpos);

            pool.install(|| rayon::join(
                || Model::parallel_copy(o1, n1, numthreads * 2, maxthreads, pool),
                || Model::parallel_copy(o2, n2, numthreads * 2, maxthreads, pool)
            ));
        } else {
            valold.copy_from_slice(valnew);
        }
    }



    fn store_data(old_temp:&mut [f64], new_temp:&[f64], datastorer:&mut [Vec<f64>], min_temp_diff: &f64, global_time: &f64, numthreads:usize, maxthreads: usize, pool: &ThreadPool){
        if numthreads <= maxthreads{
            let splitpos = old_temp.len()/2;
            let (ot1, ot2) = old_temp.split_at_mut(splitpos);
            let (nt1, nt2) = new_temp.split_at(splitpos);
            let (ds1, ds2) = datastorer.split_at_mut(splitpos);
            pool.install(||rayon::join(
                ||Model::store_data(ot1, nt1, ds1, min_temp_diff, global_time, numthreads*2, maxthreads, pool),
                ||Model::store_data(ot2, nt2, ds2, min_temp_diff, global_time, numthreads*2, maxthreads, pool)
            ));
        }
        else{
            for i in 0..old_temp.len(){
                if (new_temp[i]-old_temp[i]).abs() > *min_temp_diff {
                    datastorer[i].push(*global_time);
                    datastorer[i].push(new_temp[i]);
                    old_temp[i] = new_temp[i];
                }
            }
        }
    }

    fn store_data_node(old_temp:&[f64], new_temp:&[f64], datastorernode:&mut [Vec<f64>], nd_to_elem:&[Vec<usize>], nodetemp: &mut [f64], nodetempold: &mut[f64],min_temp_diff: &f64, global_time: &f64, numthreads:usize, maxthreads: usize, pool: &ThreadPool){
        if numthreads <= maxthreads{
            let splitpos = datastorernode.len()/2;
            //println!("{}",datastorernode.len());
            let (ds1, ds2) = datastorernode.split_at_mut(splitpos);
            let (nte1, nte2) = nd_to_elem.split_at(splitpos);
            let (nt1, nt2) = nodetemp.split_at_mut(splitpos);
            let (ntold1, ntold2) = nodetempold.split_at_mut(splitpos);
            pool.install(||rayon::join(
                ||Model::store_data_node(old_temp, new_temp, ds1, nte1 , nt1, ntold1,min_temp_diff, global_time, numthreads*2, maxthreads, pool),
                ||Model::store_data_node(old_temp, new_temp, ds2, nte2, nt2, ntold2, min_temp_diff, global_time, numthreads*2, maxthreads, pool)
            ));
        }
        else{
            for i in 0..datastorernode.len(){
                let mut updated = false;
                let elems_to_check = &nd_to_elem[i];
                // for j in 0..elems_to_check.len(){
                //    if elems_to_check[j]< old_temp.len() && (new_temp[elems_to_check[j]] - old_temp[elems_to_check[j]]).abs() > * min_temp_diff{
                //       updated = true;
                //       break;
                //   }
                //  }
                // if updated{
                let mut temp = 0.0;
                let mut numelems = 0;
                for j in 0..elems_to_check.len(){
                    if elems_to_check[j]<new_temp.len() {
                        temp = temp + new_temp[elems_to_check[j]];
                        numelems += 1;
                    }
                }
                temp = temp/numelems as f64;
                nodetemp[i] = temp;
                if (temp - nodetempold[i]).abs() > *min_temp_diff {
                    datastorernode[i].push(*global_time);
                    datastorernode[i].push(temp);
                    nodetempold[i] = temp;
                }
                // }
            }
        }
    }

    /// updates model based on input. Needs to be placed in Model implementation later.
    pub fn update_model(mu: &mut ModelUpdater, mdl: &mut Model, bw: &mut BufWriter<File>, areas_and_dists: [f64; 5], time_step: f64,
                        maxthreads: usize, input_data: ([f64; 3], f64, f64, f64, [f64; 2]),
                        temp_diff: f64, tempstorer: &mut Vec<Vec<f64>>,
                        tempnodestorer: &mut Vec<Vec<f64>>,
                        nd_to_elems_upd: &[Vec<usize>],

                        cooldown_period: f64,
                        pool: &ThreadPool) {
        //let max_cpus = num_cpus::get_physical();

        //let pool = rayon::ThreadPoolBuilder::new().num_threads(maxthreads).build().unwrap();
        //lists for interpolation
        let mut temporary_templist = Vec::with_capacity(mu.activation_times.len());


        let mut cell_indices = Vec::with_capacity(mu.activation_times.len());
        //println!("activation times len {}" ,mu.activation_times.len());
        let mut next_cell_info = mu.get_next_cell_info();
        let mut next_cell_info2 = mu.get_next_cell_info();
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

        println!("kx ky kz  {} {} {}", kx, ky, kz);;
        println!("t_env {}", t_env);

        mdl.addCell(next_cell_info.1, sp_heat_cap, density, kx, ky, kz, init_temp,
                    next_cell_info.3, next_cell_info.2, next_cell_info.4, &mut temporary_templist,
                    next_cell_info.5, next_cell_info.6, &mut cell_indices,  tempnodestorer);
        //global_time = global_time + time;
        global_time = global_time + time;
        if time>1e-4 {
            mdl.run_model(time, dt, conv_coeff, t_env, bw, &mut global_time, &pool, areas_and_dists, &mut temporary_templist,
                          maxthreads, &cell_indices);
        }
        Model::store_data_node(&mdl.oldtemplist, &temporary_templist, tempnodestorer,nd_to_elems_upd, &mut mdl.nodetemp, &mut mdl.nodetemp_old, &temp_diff, &global_time,1, maxthreads,pool);
        Model::store_data(&mut mdl.oldtemplist, &temporary_templist, tempstorer,&temp_diff, &global_time,1, maxthreads,pool);


        //println!("run model");
        let endval = mu.activation_times.len();
        for i in 1..endval - 1 {
            next_cell_info = next_cell_info2;
            next_cell_info2 = mu.get_next_cell_info();

            let time = next_cell_info2.0 - next_cell_info.0;
            // global_time = global_time + time;
            if i % 1000 == 0 { println!("{} out of {}, time {}", i, endval, global_time);}
            //println!("layer {}",next_cell_info.6);
            mdl.addCell(next_cell_info.1, sp_heat_cap, density, kx, ky, kz, init_temp,
                        next_cell_info.3, next_cell_info.2, next_cell_info.4, &mut temporary_templist,
                        next_cell_info.5, next_cell_info.6, &mut cell_indices, tempnodestorer);
            //println!("second addcell");
            //println!("second run_model");
            //global_time = global_time + time;
            //println!("time {}", time);
            let timespace = dt;
            // if time>10.0 {println!("time > 10, global time {}", global_time);}
            // else {println!("time is less than 10, global time {}", global_time);}
            if time > 1e-4 {

                let nsteps = (time/timespace) as usize;
                for i in 0..nsteps {
                    global_time = global_time + timespace;
                    mdl.run_model(timespace, dt, conv_coeff, t_env, bw, &mut global_time, &pool, areas_and_dists, &mut temporary_templist,
                                  maxthreads, &cell_indices);
                    Model::store_data_node(&mdl.oldtemplist, &temporary_templist, tempnodestorer,nd_to_elems_upd,  &mut mdl.nodetemp, &mut mdl.nodetemp_old,&temp_diff, &global_time,1, maxthreads,pool);
                    // println!("global time {}", global_time);
                    // println!("time timespace nsteps {} {} {}", time, timespace, nsteps);
                    Model::store_data(&mut mdl.oldtemplist, &temporary_templist, tempstorer,&temp_diff, &global_time,1, maxthreads,pool);
                }
                let remaining_time = time - timespace * (nsteps as f64);
                global_time = global_time + remaining_time;
                mdl.run_model(timespace, dt, conv_coeff, t_env, bw, &mut global_time, &pool, areas_and_dists, &mut temporary_templist,
                              maxthreads, &cell_indices);
                Model::store_data_node(&mdl.oldtemplist, &temporary_templist, tempnodestorer,nd_to_elems_upd,  &mut mdl.nodetemp, &mut mdl.nodetemp_old,&temp_diff, &global_time,1, maxthreads,pool);
                Model::store_data(&mut mdl.oldtemplist, &temporary_templist, tempstorer,&temp_diff, &global_time,1, maxthreads,pool);
            }

            //Model::store_data_node(&mdl.oldtemplist, &temporary_templist, tempnodestorer,nd_to_elems_upd,  &mut mdl.nodetemp, &mut mdl.nodetemp_old,&temp_diff, &global_time,1, maxthreads,pool);
            //Model::store_data(&mut mdl.oldtemplist, &temporary_templist, tempstorer,&temp_diff, &global_time,1, maxthreads,pool);
            //println!("second run model");
        }
        next_cell_info = next_cell_info2;

        mdl.addCell(next_cell_info.1, sp_heat_cap, density, kx, ky, kz, init_temp,
                    next_cell_info.3, next_cell_info.2, next_cell_info.4, &mut temporary_templist,

                    next_cell_info.5, next_cell_info.6, &mut cell_indices, tempnodestorer);

        //global_time = global_time + time;
        //let timedt = 2.0;
        for i in 0..(cooldown_period / dt) as usize {
            global_time = global_time + dt;
            mdl.run_model(time, dt, conv_coeff, t_env, bw, &mut global_time, &pool, areas_and_dists, &mut temporary_templist,
                          maxthreads, &cell_indices);

            Model::store_data_node(&mdl.oldtemplist, &temporary_templist, tempnodestorer,nd_to_elems_upd,  &mut mdl.nodetemp, &mut mdl.nodetemp_old,&temp_diff, &global_time,1, maxthreads,pool);
            Model::store_data(&mut mdl.oldtemplist, &temporary_templist, tempstorer,&temp_diff, &global_time,1, maxthreads,pool);
        }

        for i in 0..tempstorer.len(){
            tempstorer[i].push(global_time);
            tempstorer[i].push(temporary_templist[i]);
        }
        println!("{} out of {}, time {}", endval, endval, global_time);
    }
}
