use crate::primitives::{Node, Point};
use rayon::ThreadPool;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, BufReader, BufRead};
use std::collections::HashMap;
use std::{fmt, fs};
use std::io::Write;
use std::str::FromStr;

/// ModelGenerator generates the domain with nodes, connectivity, and neighbor information
#[derive(Clone, Debug)]
pub struct ModelGenerator{
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
    zmin: f64,
    zmax: f64,
    x_divisions: usize,
    y_divisions: usize,
    z_divisions: usize,
    pub(crate) nodelist: Vec<Node>,
    pub(crate) celllist: Vec<[usize;8]>,
    pub(crate) neighborlist:Vec<[isize;6]>,
    pub(crate) has_neighbor:Vec<[bool;6]>,
    element_width:f64,
    centers:Vec<Point>,
    z_divs_per_layer: usize,
}

/// methods for ModelGenerator
impl ModelGenerator{
    pub fn new(xmin:f64, xmax:f64, ymin:f64, ymax:f64, zmin:f64,zmax:f64, x_divisions:usize, y_divisions:usize, z_divisions:usize, init_temp:f64, element_width:f64, z_divs_per_layer: usize)->ModelGenerator{
        let mut mdlgen = ModelGenerator::new_internal(xmin, xmax,ymin,ymax,zmin,zmax,x_divisions,y_divisions,z_divisions,init_temp, element_width, z_divs_per_layer);
        mdlgen.calc_centers();
        return mdlgen
    }

    pub fn save_as_file(&self, filename:&str){
        let mut base_filename = filename.to_string();
        let mut fn1 = base_filename.clone(); fn1.push('1');
        let mut fn2 = base_filename.clone(); fn2.push('2');
        let mut fn3 = base_filename.clone(); fn3.push('3');
        let mut fn4 = base_filename.clone(); fn4.push('4');
        let mut fn5 = base_filename.clone(); fn5.push('5');
        let mut fn6 = base_filename.clone(); fn6.push('5');

        let mut file1 = File::create(fn1).expect("can't create the file");
        let mut file2 = File::create(fn2).expect("can't create the file");
        let mut file3 = File::create(fn3).expect("can't create the file");
        let mut file4 = File::create(fn4).expect("can't create the file");
        let mut file5 = File::create(fn5).expect("can't create the file");
        let mut file6 = File::create(fn6).expect("can't create the file");

        let mut bf1 = BufWriter::with_capacity(100000,file1);
        let mut bf2 = BufWriter::with_capacity(100000,file2);
        let mut bf3 = BufWriter::with_capacity(100000,file3);
        let mut bf4 = BufWriter::with_capacity(100000,file4);
        let mut bf5 = BufWriter::with_capacity(100000,file5);
        let mut bf6 = BufWriter::with_capacity(100000,file6);


        //store small stuff in file 1
        write!(bf1, "{},{},{},{},{},{},{},{},{},{},{}",self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax,self.x_divisions, self.y_divisions, self.z_divisions, self.element_width, self.z_divs_per_layer);
        //store nodelist in bf2
        for i in &self.nodelist{
            writeln!(bf2, "{},{},{},{}",i.index,i.pt.x, i.pt.y, i.pt.z);
        }
        //store celllist in bf3
        for i in &self.celllist{
            writeln!(bf3, "{},{},{},{},{},{},{},{}",i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7]);
        }
        //store has_neighbors in bf4
        for i in &self.neighborlist{
            writeln!(bf4, "{},{},{},{},{},{}",i[0],i[1],i[2],i[3],i[4],i[5]);
        }
        //store neighbors in bf5
        for i in &self.has_neighbor{
            writeln!(bf5, "{},{},{},{},{},{}",i[0],i[1],i[2],i[3],i[4],i[5]);
        }
        //store centers in bf6
        for i in &self.centers{
            writeln!(bf6,"{},{},{}",i.x, i.y, i.z);
        }
    }

    pub fn from_files(filename: &str)->ModelGenerator{
        let mut base_filename = filename.to_string();
        let mut fn1 = base_filename.clone(); fn1.push('1');
        let mut fn2 = base_filename.clone(); fn2.push('2');
        let mut fn3 = base_filename.clone(); fn3.push('3');
        let mut fn4 = base_filename.clone(); fn4.push('4');
        let mut fn5 = base_filename.clone(); fn5.push('5');
        let mut fn6 = base_filename.clone(); fn6.push('5');


        let mut file1 = File::open(fn1).expect("can't open the file");
        let mut file2 = File::open(fn2).expect("can't open the file");
        let mut file3 = File::open(fn3).expect("can't open the file");
        let mut file4 = File::open(fn4).expect("can't open the file");
        let mut file5 = File::open(fn5).expect("can't open the file");
        let mut file6 = File::open(fn6).expect("can't open the file");

        let mut bf1 = BufReader::with_capacity(100000,file1);
        let mut bf2 = BufReader::with_capacity(100000,file2);
        let mut bf3 = BufReader::with_capacity(100000,file3);
        let mut bf4 = BufReader::with_capacity(100000,file4);
        let mut bf5 = BufReader::with_capacity(100000,file5);
        let mut bf6 = BufReader::with_capacity(100000,file6);

        let mut temp = String::new();
        bf1.read_line(&mut temp).expect("cant read line");
        let mut all_dat = temp.split(",");
        let xmin = f64::from_str(all_dat.next().expect("no next")).expect("can't parse");
        let xmax = f64::from_str(all_dat.next().expect("no next")).expect("can't parse");
        let ymin = f64::from_str(all_dat.next().expect("no next")).expect("can't parse");
        let ymax = f64::from_str(all_dat.next().expect("no next")).expect("can't parse");
        let zmin = f64::from_str(all_dat.next().expect("no next")).expect("can't parse");
        let zmax = f64::from_str(all_dat.next().expect("no next")).expect("can't parse");
        let x_divisions = usize::from_str(all_dat.next().expect("no next")).expect("can't parse");
        let y_divisions = usize::from_str(all_dat.next().expect("no next")).expect("can't parse");
        let z_divisions = usize::from_str(all_dat.next().expect("no next")).expect("can't parse");
        let element_width = f64::from_str(all_dat.next().expect("no next")).expect("can't parse");
        let z_divs_per_layer = usize::from_str(all_dat.next().expect("no next")).expect("can't parse");

        let mut nodelist = Vec::with_capacity(10000);
        for i in bf2.lines(){
            let mut tmpi = i.expect("can't read string");
            let mut tmp = tmpi.split(",");
            let index = usize::from_str(tmp.next().expect("no next")).expect("can't parse");
            let x = f64::from_str(tmp.next().expect("no next")).expect("cant parse");
            let y = f64::from_str(tmp.next().expect("no next")).expect("cant parse");
            let z = f64::from_str(tmp.next().expect("no next")).expect("cant parse");
            nodelist.push(Node{ pt: Point { x, y, z }, index });
        }

        let mut celllist = Vec::with_capacity(10000);
        for i in bf3.lines(){
            let mut tmpi = i.expect("can't read string");
            let mut tmp = tmpi.split(",");
            let mut elems = [0;8];
            for j in 0..8 {
                elems[j] = usize::from_str(tmp.next().expect("no next")).expect("can't parse");
            }
            celllist.push(elems);
        }

        let mut neighborlist = Vec::with_capacity(10000);
        for i in bf4.lines(){
            let mut tmpi = i.expect("can't read string");
            let mut tmp = tmpi.split(",");
            let mut elems = [0;6];
            for j in 0..6 {
                elems[j] = isize::from_str(tmp.next().expect("no next")).expect("can't parse");
            }
            neighborlist.push(elems);
        }

        let mut has_neighbor = Vec::with_capacity(10000);
        for i in bf5.lines(){
            let mut tmpi = i.expect("can't read string");
            let mut tmp = tmpi.split(",");
            let mut elems = [false;6];
            for j in 0..8 {
                elems[j] = bool::from_str(tmp.next().expect("no next")).expect("can't parse");
            }
            has_neighbor.push(elems);
        }

        let mut centers = Vec::with_capacity(10000);
        for i in bf6.lines(){
            let mut tmpi = i.expect("can't read string");
            let mut tmp = tmpi.split(",");
            //let mut elems = [false;6];
            let x = f64::from_str(tmp.next().expect("no next")).expect("cant parse to f64");
            let y = f64::from_str(tmp.next().expect("no next")).expect("cant parse to f64");
            let z = f64::from_str(tmp.next().expect("no next")).expect("cant parse to f64");
            centers.push(Point{x,y,z});
        }

        return ModelGenerator{
            xmin,
            xmax,
            ymin,
            ymax,
            zmin,
            zmax,
            x_divisions,
            y_divisions,
            z_divisions,
            nodelist,
            celllist,
            neighborlist,
            has_neighbor,
            element_width,
            centers,
            z_divs_per_layer
        }

    }

    fn new_internal(xmin:f64, xmax:f64, ymin:f64, ymax:f64, zmin:f64,zmax:f64, x_divisions:usize, y_divisions:usize, z_divisions:usize, init_temp:f64, element_width: f64, z_divs_per_layer: usize)->ModelGenerator{
        //let z_divisions = z_divisions+1;
        let nodelist = ModelGenerator::generate_nodelist(xmin,xmax,ymin,ymax,zmin,zmax,x_divisions,y_divisions,z_divisions);
        let celllist = ModelGenerator::generate_celllist(x_divisions,y_divisions,z_divisions);
        let neighborlist = ModelGenerator::generate_neighborslist(x_divisions,y_divisions,z_divisions);
        //let init_templist = ModelGenerator::generate_initial_temperaturelist(x_divisions,y_divisions,z_divisions,init_temp);
        let centers = Vec::with_capacity((x_divisions-1)*(y_divisions-1));

        return ModelGenerator{
            xmin,
            xmax,
            ymin,
            ymax,
            zmin,
            zmax,
            x_divisions,
            y_divisions,
            z_divisions,
            nodelist,
            celllist,
            neighborlist:neighborlist.1,
            has_neighbor:neighborlist.0,
            centers,
            element_width,
            z_divs_per_layer,
        }
    }

    ///generate list of nodes sorted by z, then by y, then by x
    pub fn generate_nodelist(xmin:f64, xmax:f64, ymin:f64,ymax:f64,zmin:f64,zmax:f64,x_divisions:usize,y_divisions:usize,z_divisions:usize)-> Vec<Node>{
        let mut pointlist = Vec::with_capacity(x_divisions * y_divisions * z_divisions);
        let mut xvals = Vec::with_capacity(x_divisions);
        let mut yvals = Vec::with_capacity(y_divisions);
        let mut zvals = Vec::with_capacity(z_divisions);

        for i in 0 .. x_divisions+1{
            xvals.push( xmin+(xmax-xmin)/((x_divisions-1) as f64) *i as f64);
        }

        for i in 0..y_divisions+1{
            yvals.push(ymin+(ymax-ymin)/((y_divisions-1) as f64) * i as f64);
        }

        for i in 0..z_divisions+1{
            zvals.push(zmin+(zmax-zmin)/((z_divisions-1) as f64) * i as f64);
        }

        for i in 0..x_divisions{
            for j in 0..y_divisions{
                for k in 0..z_divisions{
                    let index = ((i)*y_divisions+j)*z_divisions+k;
                    let pt = Point{
                        x:xvals[i],
                        y:yvals[j],
                        z:zvals[k],
                    };
                    pointlist.push(pt);
                }
            }
        }
        pointlist.sort_by(|a,b| (a.x).partial_cmp(&(b.x)).unwrap());
        pointlist.sort_by(|a,b| (a.y).partial_cmp(&(b.y)).unwrap());
        pointlist.sort_by(|a,b| (a.z).partial_cmp(&(b.z)).unwrap());

        let mut nodelist = Vec::with_capacity(pointlist.len());
        for i in 0..pointlist.len().clone(){
            let index = i;
            let nd = Node{index:index, pt:pointlist[i]};
            nodelist.push(nd);
        }

        return nodelist;
    }
    /// generate list of cells with connectivity based on nodelist. Nodes for cells are in order
    /// [0,1,3,2,4,5,7,6], 0,1,3,2 lower nodes continuous in anti-clockwise direction, 4,5,7,6 upper
    /// nodes continuous in anti-clockwise direction
    pub fn generate_celllist(x_divisions:usize, y_divisions:usize, z_divisions:usize)->Vec<[usize;8]>{
        let mut celllist=Vec::with_capacity(x_divisions*y_divisions*z_divisions);
        let xincr = x_divisions;
        let yincr = y_divisions;
        let xyincr = xincr*yincr;
        for i in 0..x_divisions-1{
            celllist.push([i,i+1,i+xincr,i+1+xincr,
                i+xyincr,i+1+xyincr, i+xincr+xyincr, i+1+xincr+xyincr]);
        }

        let lencelllist = celllist.len();
        for j in 1..y_divisions-1{
            for i in 0..lencelllist{
                let mut newcell =[0,0,0,0,0,0,0,0];
                for k in 0..8{
                    newcell[k] = celllist[i][k]+xincr*j;
                }
                celllist.push(newcell);
            }
        }
        let lencelllist = celllist.len();
        for j in 1..z_divisions-1{
            for i in 0..lencelllist{
                let mut newcell =[0,0,0,0,0,0,0,0];
                for k in 0..8{
                    newcell[k] = celllist[i][k]+xyincr*j;
                }
                celllist.push(newcell);
            }
        }
        return celllist;
    }

    /// returns a list of index of neighbors. Neighbors are [front, back, left, right, up, down] or
    /// [+x, -x, +y, -y, +z, -z]. If there is no neighbor, the boolean vector gives false, and the
    /// index vector gives negative value
    pub fn generate_neighborslist(x_divisions:usize, y_divisions:usize, z_divisions:usize)->(Vec<[bool;6]>,Vec<[isize;6]>){
        let xincr = x_divisions-1;
        let yincr = y_divisions-1;
        let xyincr = xincr * yincr;
        let num_cells = (x_divisions-1) * (y_divisions-1) * (z_divisions-1);
        let mut out_neighbors = Vec::with_capacity(num_cells);
        let mut has_neighbors= Vec::with_capacity(num_cells);
        //let mut out_neighbors_temp = Vec::with_capacity(num_cells);
        let xincrt = (xincr) as isize;
        let yincrt = (yincr) as isize;
        let xyincrt = xincrt * yincrt as isize;
        for i in 0..num_cells{
            let t = i as isize;
            out_neighbors.push([(t+1), (t-1) ,(t+xincrt) ,(t -xincrt), (t+xyincrt), (t -xyincrt)]);
        }

        for i in 0..num_cells{
            let mut bool_neighbor = [true,true,true,true,true,true];
            if ((i +1) % xincr == 0) {bool_neighbor[0] = false} // has no neighbor in front (+x dir)
            if (i % xincr == 0) {bool_neighbor[1] = false} // has no neighbor in back (-x dir)
            if (xyincr - i % xyincr <= xincr ) {bool_neighbor[2] = false} // has no neighbor in left (+y dir)
            if (i % xyincr < xincr) {bool_neighbor[3] = false} // has no neighbor in right (-y dir)
            if (i >= xyincr*(z_divisions-2)) {bool_neighbor[4] = false} // has no neighbor on top (+z dir)
            if (i < xyincr) {bool_neighbor[5] = false} // has no neighbor on bottom (-z dir)
            has_neighbors.push(bool_neighbor);
        }
        return (has_neighbors, out_neighbors)
    }

    pub fn generate_initial_temperaturelist(x_divisions:usize,y_divisions:usize, z_divisions:usize,extrusion_temp:f64)->Vec<f64>{
        let num_cells = (x_divisions-1) * (y_divisions-1) * (z_divisions-1);
        let mut out_temp = Vec::with_capacity(num_cells);
        for i in 0..num_cells{
            out_temp.push(extrusion_temp);
        }
        return out_temp;
    }
    /// groups segments by layers so that each layer can be independently processed from another
    pub fn group_segments_by_layers(&self, segments:Vec<[Point;2]>, is_extrusion_on:Vec<bool>, move_speed:Vec<f64>)->(Vec<Vec<[Point;2]>>,Vec<Vec<bool>>, Vec<Vec<f64>>){
        let numlayers = (self.z_divisions-1)/self.z_divs_per_layer;
        //println!("z_divs {}, z_divs_per_layer {}, numlayers {}", self.z_divisions, self.z_divs_per_layer, numlayers);
        let mut segment_grp = Vec::with_capacity(numlayers);
        let mut is_extrusion_on_grp = Vec::with_capacity(numlayers);
        let mut move_speed_grp = Vec::with_capacity(numlayers);

        //println!("numlayers = {}", numlayers);
        for i in 0..numlayers{
            segment_grp.push(Vec::with_capacity(1000));
            is_extrusion_on_grp.push(Vec::with_capacity(1000));
            move_speed_grp.push(Vec::with_capacity(1000));
        }

        let mut oldz = segments[0][0].z;
        let mut layer_counter = 0;
        //let mut count_z_div_in_layer = 0;
        for i in 0 .. segments.len(){
            let mut newz = segments[i][0].z;
            let layer_change = (oldz-newz).abs()>1e-8;
            if layer_change {
                oldz = newz;
                newz = segments[i][0].z;
                layer_counter = layer_counter+1;
            }
            //let zz = segments[i].clone();
            //           println!("layer_counter {}, num_layers {}",layer_counter, numlayers );
            segment_grp[layer_counter].push(segments[i].clone());
            is_extrusion_on_grp[layer_counter].push(is_extrusion_on[i].clone());
            move_speed_grp[layer_counter].push(move_speed[i]);
        }

        return (segment_grp, is_extrusion_on_grp, move_speed_grp);
    }
    pub fn group_cells_by_layers(&self,beadwidth:f64 ,beadheight:f64)->Vec<Vec<usize>>{
        let num_vertical_elems = self.z_divisions-1;
        let num_layers = num_vertical_elems/self.z_divs_per_layer;

        let num_cells_in_layer = (self.x_divisions-1)*(self.y_divisions-1)*self.z_divs_per_layer;
        //println!("num layers {} num vertical elems {}, num cells in layers {}", num_layers, num_vertical_elems, num_cells_in_layer);
        // let layer_height = (self.zmax-self.zmin)/(self.z_divisions-1)as f64;
        let mut cells_grp = Vec::with_capacity(num_layers);
        for i in 0..num_layers{
            cells_grp.push(Vec::with_capacity(num_cells_in_layer));
        }
        for i in 0 .. num_layers{
            for j in 0..num_cells_in_layer{
                cells_grp[i].push(i*num_cells_in_layer+j);
            }
        }
        // println!("l1 {}, l2 {}, l3 {},l4 {}", cells_grp[0][0],cells_grp[1][0],cells_grp[2][0], cells_grp[4][0]);
        return cells_grp;
    }

    ///generate activation times for elements in a single layer
    pub fn generate_activation_times_for_layer(&self, segments:&Vec<[Point;2]>, is_extrusion_on:&Vec<bool>,
                                               move_speed:&Vec<f64>,bead_width:f64, bead_height:f64,
                                               celllist_layer:&Vec<usize>, layer_no:usize, num_cells_in_layer:usize,
                                               is_activated:&mut[bool], activation_time:&mut[f64], orientation_x:&mut[f64],
                                               orientation_y:&mut[f64], layer_no_out:&mut[usize],cell_no:&mut[usize], layer_end_time:&mut f64)
    {
        //println!("layer no {}", layer_no);
        //println!("layer no {}", layer_no);
        let width = bead_width/2.0;
        //let num_layers = self.z_divisions-1;
        //println!("num cells in layer {}", num_cells_in_layer);
        let mut is_cell_activated = Vec::with_capacity(num_cells_in_layer);
        //println!("num cells in layer {}", num_cells_in_layer);
        //let num_cells_in_layer = (self.x_divisions-1)*(self.y_divisions-1)*self.z_divs_per_layer;
        /*let mut centers = Vec::with_capacity(num_cells_in_layer);
        for i in 0..num_cells_in_layer{
            centers.push(self.calc_center(self.celllist[i]));
        }*/

        for i in 0..num_cells_in_layer{
            is_cell_activated.push(false);
        }
        let mut time = 0.0 as f64;
        for i in 0..segments.len() {
            let current_segment = segments[i];
            let dist = current_segment[0].dist_to(&current_segment[1]);
            let dt = dist / move_speed[i];
            let xdir = current_segment[1].x - current_segment[0].x;
            let ydir = current_segment[1].y - current_segment[0].y;

            let length = (xdir * xdir + ydir * ydir).sqrt();
            let orient = [xdir/length, ydir/length];
            if is_extrusion_on[i] {
                let divs = (length/width) as usize + 2;
                let time_and_pos_vec = ModelGenerator::get_pos_and_time(time, current_segment[0], current_segment[1], divs, move_speed[i]);
                // check only cells within certain width from the bead being deposited. This is a filter to speed up the actual checking
                let mut cells_to_check = Vec::with_capacity(100000);
                let y1 = current_segment[0].y; let y2 = current_segment[1].y;
                let x1 = current_segment[0].x; let x2 = current_segment[1].x;
                let m = (y2 - y1) / (x2 - x1); let c = y1 - m * x1;

                for ii in 0..num_cells_in_layer{
                    let center = self.centers[ii];
                    let centery_predicted = m * center.x + c;
                    let err;
                    if (x2-x1).abs()>1e-6 { // if not equal x1, x2
                        err = (center.y-centery_predicted).abs();
                    }
                    else{
                        err = (center.x-x2).abs();
                    }
                    if (err < 12.0 * bead_width) {cells_to_check.push(ii);} // need to check why 12 is the magic number??
                }
                for j in time_and_pos_vec {
                    for k in 0..cells_to_check.len() {
                        let current_cell = cells_to_check[k] ;
                        if !is_cell_activated[current_cell] {

                            //  if self.segment_intersects_cell(&self.celllist[current_cell],[current_pt,next_pt],0.8,0.4)
                            let center = self.calc_center(self.celllist[current_cell]);
                            //let center2 = j.1;
//                            let width = bead_width/2.0;
                            let width_small = 1.0 * width; //reduce here to to check within smaller area
                            let p1 = Point { x: j.1.x - width_small, y: j.1.y - width_small, z: j.1.z };
                            let p2 = Point { x: j.1.x + width_small, y: j.1.y - width_small, z: j.1.z };
                            let p3 = Point { x: j.1.x + width_small, y: j.1.y + width_small, z: j.1.z };
                            let p4 = Point { x: j.1.x - width_small, y: j.1.y + width_small, z: j.1.z };


                            if ModelGenerator::center_lies_inside_area([p1, p2, p3, p4], center) & !is_cell_activated[current_cell]
                            {
                                is_cell_activated[current_cell] = true;

                                let this_cell_in_slice = current_cell;
                                activation_time[this_cell_in_slice] = j.0;
                                cell_no[this_cell_in_slice] = current_cell+num_cells_in_layer*layer_no;
                                layer_no_out[this_cell_in_slice] = layer_no;
                                orientation_x[this_cell_in_slice] = orient[0];
                                orientation_y[this_cell_in_slice] = orient[1];
                                is_activated[this_cell_in_slice] = true;

                            }
                        }
                    }

                }

            }
            time = time + dt;
        }
        *layer_end_time = time;
    }

    fn calc_centers(&mut self){
        for i in 0..(self.x_divisions-1)*(self.y_divisions-1)*self.z_divs_per_layer{
            self.centers.push(self.calc_center(self.celllist[i]));
        }
    }

    pub fn generate_activation_times_all_layers(&self,segments:Vec<[Point;2]>, is_extrusion_on:Vec<bool>,
                                                move_speed:Vec<f64>,bead_width:f64, bead_height:f64,
                                                pool:&ThreadPool, max_threads:usize) -> Vec<(f64,usize,[f64;2],usize)>
    {

        //self.calc_centers();
        let grps_segs = self.group_segments_by_layers(segments,is_extrusion_on,move_speed);
        let grps_cells = self.group_cells_by_layers(bead_width, bead_height);
        //&self, grps_segs:&(Vec<Vec<[Point;2]>>,Vec<Vec<bool>>, Vec<Vec<f64>>), grps_cells:&Vec<Vec<usize>>,
        //     beadwidth:f64, beadheight:f64
        let num_layers = (self.z_divisions - 1)/self.z_divs_per_layer;
        let num_cells_in_layer = (self.x_divisions - 1) * (self.y_divisions - 1) * self.z_divs_per_layer;
        //println!("num cells in layer {}", num_cells_in_layer);

        let layer_numbers:Vec<usize> = (0..num_layers).map(|i|i).collect();
        // let (ln0,ln1) = layer_numbers.split_at(num_layers/2);

        //
        let segments_grp = &grps_segs.0;

        let extrusion_grp = &grps_segs.1;
        let movespeed_grp = &grps_segs.2;
        //let (sg0,sg1) = segments_grp.split_at(num_layers/2);
        // let(ext0,ext1) = extrusion_grp.split_at(num_layers/2);
        // let(mvs0, mvs1) = movespeed_grp.split_at(num_layers/2);



        let (grpc0,grpc1) = grps_cells.split_at(num_layers/2);

        let mut layer_end_time_grp = vec![0.0 as f64; num_layers];
        let (leg0, leg1) = layer_end_time_grp.split_at_mut(num_layers/2);

        let mut is_activated_grp = Vec::with_capacity(num_layers);
        let mut activation_time_grp = Vec::with_capacity(num_layers);
        let mut orientation_x_grp = Vec::with_capacity(num_layers);
        let mut orientation_y_grp = Vec::with_capacity(num_layers);
        let mut no_layer_out_grp = Vec::with_capacity(num_layers);
        let mut cell_no_grp = Vec::with_capacity(num_layers);

        for i in 0..num_layers {
            is_activated_grp.push(vec![false; num_cells_in_layer]);
            activation_time_grp.push(vec![-1.0 as f64; num_cells_in_layer]);
            orientation_x_grp.push(vec![f64::NEG_INFINITY; num_cells_in_layer]);
            orientation_y_grp.push(vec![f64::NEG_INFINITY; num_cells_in_layer]);
            no_layer_out_grp.push(vec![0 as usize; num_cells_in_layer]);
            cell_no_grp.push(vec![0 as usize; num_cells_in_layer]);
        }

        self.calc_act_tim_split(segments_grp, extrusion_grp, movespeed_grp, &grps_cells, layer_end_time_grp.as_mut_slice(),
                                is_activated_grp.as_mut_slice(),activation_time_grp.as_mut_slice(),
                                orientation_x_grp.as_mut_slice(), orientation_y_grp.as_mut_slice(),
                                no_layer_out_grp.as_mut_slice(), cell_no_grp.as_mut_slice(), 1, max_threads,
                                layer_numbers.as_slice(),num_layers, bead_width, bead_height,pool,
                                num_cells_in_layer );

        let mut activation_times = Vec::with_capacity(self.celllist.len() / 2);
        let mut previous_layer_time = 0.0 as f64;
        for i in 0..num_layers {
            for j in 0..num_cells_in_layer {
                if is_activated_grp[i][j] {
                    activation_times.push((activation_time_grp[i][j]+previous_layer_time,
                                           cell_no_grp[i][j],[orientation_x_grp[i][j], orientation_y_grp[i][j]], i));
                }
            }
            previous_layer_time = previous_layer_time+layer_end_time_grp[i];
        }
        activation_times.sort_by(|a,b|PartialOrd::partial_cmp(&a.0,&b.0).expect("cant compare these values, either NaN or Infinity"));
        return activation_times;
    }

    fn calc_act_tim_split(&self, segments_grp:&[Vec<[Point;2]>], extrusion_grp: &[Vec<bool>], movespeed_grp:&[Vec<f64>],
                          grps_cells:&[Vec<usize>],layer_time_end_grp:&mut [f64], is_activated_grp:&mut [Vec<bool>],
                          activation_time_grp:&mut [Vec<f64>], orientation_x_grp:&mut [Vec<f64>], orientation_y_grp:&mut [Vec<f64>],
                          no_layer_out_grp:&mut [Vec<usize>], cell_no_grp: &mut [Vec<usize>], numthreads: usize, maxthreads:usize,
                          layer_numbers:&[usize],num_layers:usize, beadwidth:f64, beadheight:f64,pool:&ThreadPool, num_cells_in_layer:usize)
    {
        if numthreads < maxthreads {
            let (ln0, ln1) = layer_numbers.split_at(num_layers / 2);
            let (sg0, sg1) = segments_grp.split_at(num_layers / 2);
            let (ext0, ext1) = extrusion_grp.split_at(num_layers / 2);
            let (mvs0, mvs1) = movespeed_grp.split_at(num_layers / 2);
            let (grpc0, grpc1) = grps_cells.split_at(num_layers / 2);

            let (iag0, iag1) = is_activated_grp.split_at_mut(num_layers / 2);
            let (atg0, atg1) = activation_time_grp.split_at_mut(num_layers / 2);
            let (otx0, otx1) = orientation_x_grp.split_at_mut(num_layers / 2);
            let (oty0, oty1) = orientation_y_grp.split_at_mut(num_layers / 2);
            let (nlo0, nlo1) = no_layer_out_grp.split_at_mut(num_layers / 2);
            let (cng0, cng1) = cell_no_grp.split_at_mut(num_layers / 2);
            let (leg0, leg1) = layer_time_end_grp.split_at_mut(num_layers / 2);

            pool.install(|| rayon::join(
                || self.calc_act_tim_split(sg0, ext0, mvs0,
                                           grpc0, leg0, iag0,
                                           atg0, otx0, oty0,
                                           nlo0, cng0, numthreads * 2, maxthreads,
                                           ln0, num_layers / 2, beadwidth, beadheight, pool,
                                           num_cells_in_layer),
                || self.calc_act_tim_split(sg1, ext1, mvs1,
                                           grpc1, leg1, iag1,
                                           atg1, otx1, oty1,
                                           nlo1, cng1, numthreads * 2, maxthreads,
                                           ln1, num_layers / 2, beadwidth, beadheight, pool,
                                           num_cells_in_layer),
            ));
        } else {
            for i in 0..grps_cells.len() {
                self.generate_activation_times_for_layer(&segments_grp[i], &extrusion_grp[i], &movespeed_grp[i],
                                                         beadwidth, beadheight, &grps_cells[i], layer_numbers[i], num_cells_in_layer, &mut is_activated_grp[i],
                                                         &mut activation_time_grp[i], &mut orientation_x_grp[i], &mut orientation_y_grp[i],
                                                         &mut no_layer_out_grp[i], &mut cell_no_grp[i], &mut layer_time_end_grp[i]);
            }
        }
    }
    /// generates activation time and orientation of deposition path for each element that gets activated.
    pub fn generate_activation_times(&self,segments:Vec<[Point;2]>, is_extrusion_on:Vec<bool>, move_speed:Vec<f64>,
                                     bead_width:f64, bead_height:f64)
                                     -> Vec<(f64,usize,[f64;2],usize)>
    {
        //let mut datafile = File::create("C:\\rustFiles\\activation_point_centers.csv").expect("cant create file");
        //let mut datafile_buf = BufWriter::with_capacity(10000,datafile);
        //let mut datafile2 = File::create("C:\\rustFiles\\cell_center.csv").expect("cant create file");
        //let mut datafile_buf2 = BufWriter::with_capacity(10000,datafile2);
        //let mut datafile3 = File::create("C:\\rustFiles\\activation_pt_center_time.csv").expect("cant create file");
        //let mut datafile_buf3 = BufWriter::with_capacity(10000,datafile3);
        /* for i in &self.celllist{
             let thiscenter = self.calc_center(i.clone());
             //write!(datafile_buf2,"{},{},{}\n",thiscenter.x,thiscenter.y,thiscenter.z);
         }*/
        let width = bead_width/2.0;
        let mut activation_times = Vec::with_capacity( self.celllist.len().clone());
        //let xincr = (self.xmax-self.xmin)/self.x_divisions as f64/12.0;
        let mut time: f64 = 0.0;
        let num_cells_in_layer = (self.x_divisions-1)*(self.y_divisions-1);
        let layer_height = (self.zmax-self.zmin)/(self.z_divisions-1)as f64;
        let num_layers = self.z_divisions-1;
        let mut cell_groupings= Vec::with_capacity(num_layers);
        for i in 0..num_layers{
            let cells_in_layer = 0..num_cells_in_layer;
            let mut cells_in_layer_group = Vec::with_capacity(num_cells_in_layer);
            for j in cells_in_layer{
                cells_in_layer_group.push(j+i*num_cells_in_layer);
            }
            cell_groupings.push(cells_in_layer_group);
        }

        let mut is_cell_activated = Vec::with_capacity(num_cells_in_layer*num_layers);
        for i in 0..num_cells_in_layer*num_layers{
            is_cell_activated.push(false);
        }
// explore parallelism
// group segments by layers
        for i in 0..segments.len(){
            let current_segment = segments[i];
            let dist = current_segment[0].dist_to(&current_segment[1]);
            let dt = dist/move_speed[i];

            // divide cells into groups by layer
            if is_extrusion_on[i]{
                //write!(datafile_buf,"{},{},{},{},{},{}\n",current_segment[0].x,current_segment[0].y,current_segment[0].z,
                // current_segment[1].x,current_segment[1].y,current_segment[1].z);
                let divs = self.x_divisions*4;
                let time_and_pos_vec = ModelGenerator::get_pos_and_time(time,current_segment[0],current_segment[1],divs,move_speed[i]);
                let layer_no = ((current_segment[0].z-self.zmin)/layer_height-1.0).round() as usize;
                let minx; let maxx; let miny; let maxy;
                if current_segment[0].x>current_segment[1].x{minx = current_segment[1].x; maxx = current_segment[0].x}
                else {minx = current_segment[0].x; maxx = current_segment[1].x}
                if current_segment[0].y>current_segment[1].y{miny = current_segment[1].y; maxy = current_segment[0].y}
                else {miny = current_segment[0].y; maxy = current_segment[1].y}

                // let centers:Vec<Point> = cell_groupings[layer_no.clone()].clone().iter().map(|i|self.calc_center(self.celllist[i.clone()])).collect();

                let mut cells_to_check = Vec::with_capacity(100000);
                let y1 = current_segment[0].y; let y2 = current_segment[1].y;
                let x1 = current_segment[0].x; let x2 = current_segment[1].x;
                let m = (y2 - y1) / (x2 - x1); let c = y1 - m * x1;
                //println!("{}",i);
                for ii in cell_groupings[layer_no.clone()].clone(){
                    // println!("{},{}",self.celllist.len(),ii);
                    let center = self.calc_center(self.celllist[ii]);
                    //need to change this to match the move direction

                    let centery_predicted = m * center.x + c;
                    let err;
                    if (x2 - x1).abs() > 1e-6 {
                        err = (center.y-centery_predicted).abs();
                    }
                    else{
                        err = (center.x-x2).abs();
                    }

                    if (err < 12.0 * bead_width){ cells_to_check.push(ii);}
                    //if (center.x>minx*0.99)&(center.x<maxx*1.01)&(center.y>miny*0.99)&(center.y<maxy*1.01){ cells_to_check.push(ii);}

                }

                //let mut j0 = time_and_pos_vec[0];
                for j in time_and_pos_vec{
                    // let mut j = time_and_pos_vec[j1];
                    //write!(datafile_buf,"{},{},{},{}\n",j.1.x,j.1.y,j.1.z,j.0);

                    for k in &cells_to_check{
                        let current_cell = *k;
                        if !is_cell_activated[current_cell]{

                            //  if self.segment_intersects_cell(&self.celllist[current_cell],[current_pt,next_pt],0.8,0.4)
                            let center = self.calc_center(self.celllist[current_cell]);
                            //let center2 = j.1;

//                            let width = bead_width/2.0;
                            let p1 = Point{x:j.1.x-width,y:j.1.y-width,z:j.1.z};
                            let p2 = Point{x:j.1.x+width,y:j.1.y-width,z:j.1.z};
                            let p3 = Point{x:j.1.x+width,y:j.1.y+width,z:j.1.z};
                            let p4 = Point{x:j.1.x-width,y:j.1.y+width,z:j.1.z};
                            //let xd = center.x-center2.x;
                            //let yd = center.y-center2.y;
                            //let dist = xd*xd+yd*yd;

                            if ModelGenerator::center_lies_inside_area([p1,p2,p3,p4],center) & !is_cell_activated[current_cell]
                            {
                                is_cell_activated[current_cell] = true;

                                let xdir = current_segment[1].x-current_segment[0].x;
                                let ydir = current_segment[1].y-current_segment[0].y;

                                let divider = (xdir*xdir+ydir*ydir).sqrt();

                                activation_times.push((j.0, current_cell,[xdir/divider,ydir/divider],layer_no));
                                //write!(datafile_buf3,"{},{},{},{}\n",center.x,center.y,center.z,j.0);
                                //print!("lolol activated {},", current_cell);

                            }

                        }

                    }
                }


            }
            time = time + dt;
        }
        return activation_times;
    }

    pub fn calc_center(&self,cell_connectivity:[usize;8])->Point{
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;

        for i in 0..8{
            x = x+self.nodelist[cell_connectivity[i]].pt.x/8.0;
            y = y+self.nodelist[cell_connectivity[i]].pt.y/8.0;
            z = z+self.nodelist[cell_connectivity[i]].pt.z/8.0;
        }
        return Point{x,y,z}
    }
    pub fn get_pos_and_time(inittime:f64, pt1:Point, pt2:Point, divs:usize, speed:f64)->Vec<(f64,Point)>{
        let mut time_and_pos = Vec::with_capacity(divs);
        for i in 0..divs{
            let current_pt = Point{x:pt1.x+(pt2.x-pt1.x)/(divs-1) as f64*i as f64,
                y:pt1.y+(pt2.y-pt1.y)/(divs-1) as f64*i as f64,
                z:pt1.z+(pt2.z-pt1.z)/(divs-1) as f64*i as f64};
            let current_time = inittime+pt1.dist_to(&current_pt)/speed;
            time_and_pos.push((current_time, current_pt))
        }
        return time_and_pos
    }

    pub fn center_lies_inside_area_general(pts:[Point;4], center:Point) ->bool {
        let o1 = ModelGenerator::calc_orientation(pts[0],pts[1],center);
        let o2 = ModelGenerator::calc_orientation(pts[1],pts[2],center);
        let o3 = ModelGenerator::calc_orientation(pts[2],pts[3],center);
        let o4 = ModelGenerator::calc_orientation(pts[3],pts[0],center);

        if ((o1==o2)|(o1==0)|(o2==0)) &((o2==o3)|(o2==0)|(o3==0)) &((o3==o4)|(o3==0)|(04==0)) &((o4==o1)|(o4==0)|(o1==0)){
            return true;
        }
        else {
            return false;
        }
    }
    pub fn center_lies_inside_area(pts:[Point;4], center:Point) ->bool {
        if(center.x>pts[0].x && center.x <pts[2].x && center.y > pts[0].y && center.y<pts[2].y){
            return true;
        }
        else{
            return false;
        }
    }


    /// checks if a given cell intersects with a given line segment and if yes, the point where it
    /// intersects

    fn segment_intersects_cell(&self,cell_conenctivity:&[usize;8], segment:[Point;2])->bool{
        let intersects;
        let p1 = segment[0];
        let p2 = segment[1];

        let q1 = self.nodelist[cell_conenctivity[0]].pt;
        let q2 = self.nodelist[cell_conenctivity[1]].pt;
        let q3 = self.nodelist[cell_conenctivity[2]].pt;
        let q4 = self.nodelist[cell_conenctivity[2]].pt;


        if ModelGenerator::lines_intersect(p1,p2,q1,q2) | ModelGenerator::lines_intersect(p1,p2,q2,q3)|ModelGenerator::lines_intersect(p1,p2,q3,q4) | ModelGenerator::lines_intersect(p1,p2,q1,q4){
            intersects = true;
        }
        else{
            intersects = false;
        }

        return intersects;
    }

    pub fn lines_intersect(p1:Point,p2:Point,q1:Point,q2:Point)->bool{
        let intersects:bool;
        let o1 = ModelGenerator::calc_orientation(p1,q1,p2);
        let o2 = ModelGenerator::calc_orientation(p1,q1,q2) ;
        let o3 = ModelGenerator::calc_orientation(p2,q2,p1);
        let o4 = ModelGenerator::calc_orientation(p2,q2,q1);
        if ((o1!=o2) &  (o3!= o4)) | ((o1 == 0) & ModelGenerator::is_on_segment(p1,p2,q1)) | ((o2 ==0) &
            ModelGenerator::is_on_segment(p2,p1,q2)) | ((o3==0) & ModelGenerator::is_on_segment(p2,q1,q2)) {
            intersects = true;
        }
        else{
            intersects = false;
        }
        return intersects;
    }

    pub fn  is_on_segment (p:Point,q:Point,r:Point) -> bool {
        let is_onsegment;
        let maxx; let minx ; let maxy; let miny;
        if (p.x>r.x) {maxx = p.x; minx = r.x} else {maxx = r.x; minx = p.x}
        if (p.y>r.y) {maxy = p.y; miny = r.y} else {maxy = r.y; miny = p.y}
        if (q.x <= maxx) & (q.x >= minx) &&
            (q.y <= maxy) & (q.y >= maxx){
            is_onsegment = true;
        }
        else{
            is_onsegment = false;
        }
        return is_onsegment;
    }

    pub fn calc_orientation(p:Point,q:Point,r:Point)->u8{
        let col = ((q.y-p.y)*(r.x-q.x)-(q.x-p.x)*(r.y-q.y));
        let orient:u8;
        if col == 0.0 {
            orient = 0;
        }
        else if col > 0.0 {
            orient = 1;
        }
        else {
            orient = 2;
        }
        return orient;
    }

    pub fn write_nodes_and_elements_to_file(&self, nodefile: &str, elemfile: &str, activated_elements:Vec<usize> )->(usize,usize, Vec<[f64;3]>, Vec<[usize;8]>, HashMap<usize, usize>){
        let mut nodefl = File::create(nodefile).expect("can't create the nodefile at the specified location");
        let mut elemfl = File::create(elemfile).expect("can't create the elementfile at the specified location");
        let mut nodeflbuf = BufWriter::with_capacity(100000,nodefl);
        let mut elemflbuf = BufWriter::with_capacity(100000,elemfl);
        let mut node_table_old_to_new = HashMap::with_capacity(100000);
        let mut node_vec = Vec::with_capacity(100000);
        let mut elem_vec = Vec::with_capacity(100000);
        let mut node_counter = 0;

        for i in activated_elements{
            let elem = self.celllist[i];
            let mut elem_new:[usize;8] = [0;8];
            for j in 0..8{
                if node_table_old_to_new.contains_key(&elem[j]){
                    let new_node_number = node_table_old_to_new.get(&elem[j]).expect("cant find key");
                    elem_new[j] = *new_node_number;
                } else{
                    node_table_old_to_new.insert(elem[j].clone(),node_counter.clone());
                    elem_new[j] = node_counter.clone();
                    node_vec.push([self.nodelist[elem[j]].pt.x,self.nodelist[elem[j]].pt.y,self.nodelist[elem[j]].pt.z]);
                    node_counter += 1;
                }
            }
            elem_vec.push(elem_new);
        }
        for i in 0..node_vec.len() {
            write!(nodeflbuf, "{},{},{},{}\n",i+1,node_vec[i][0],node_vec[i][1],node_vec[i][2]);
        }

        for i in 0..elem_vec.len(){
            write!(elemflbuf,"{},{},{},{},{},{},{},{},{}\n",i+1, elem_vec[i][0]+1,elem_vec[i][1]+1,elem_vec[i][2]+1,elem_vec[i][3]+1,elem_vec[i][4]+1,elem_vec[i][5]+1,elem_vec[i][6]+1,elem_vec[i][7]+1);
        }

        return (node_vec.len(), elem_vec.len(), node_vec, elem_vec, node_table_old_to_new )


    }


    pub fn write_abaqus_input_mesh_file(&self, abaqus_meshfile: &str, activated_elements:Vec<usize>){
        let mut abaqus_input = File::create(abaqus_meshfile).expect("can't create the abaqus input file at the specified location");
        let mut abaqus_input_buffer = BufWriter::with_capacity(100000,abaqus_input);

        let mut node_table_old_to_new = HashMap::with_capacity(100000);
        let mut node_vec = Vec::with_capacity(100000);
        let mut elem_vec = Vec::with_capacity(100000);
        let mut node_counter = 0;

        for i in activated_elements{
            let elem = self.celllist[i];
            let mut elem_new:[usize;8] = [0;8];
            for j in 0..8{
                if node_table_old_to_new.contains_key(&elem[j]){
                    let new_node_number = node_table_old_to_new.get(&elem[j]).expect("cant find key");
                    elem_new[j] = *new_node_number;
                } else{
                    node_table_old_to_new.insert(elem[j].clone(),node_counter.clone());
                    elem_new[j] = node_counter.clone();
                    node_vec.push([self.nodelist[elem[j]].pt.x,self.nodelist[elem[j]].pt.y,self.nodelist[elem[j]].pt.z]);
                    node_counter += 1;
                }
            }
            elem_vec.push(elem_new);
        }

        write!(abaqus_input_buffer,
               "*Heading\n
** Job name: abaqus_meshfile Model name: Model-1\n
** Generated by: Sunil Bhandari\n
*Preprint, echo=NO, model=NO, history=NO, contact=NO\n
**\n
** PARTS\n
**\n
*Part, name=Part-1\n
*End Part\n
**\n
**\n
** ASSEMBLY\n
**\n
*Assembly, name=Assembly\n
**\n
*Instance, name=Part-1-1, part=Part-1\n
*Node\n");
        for i in 0..node_vec.len() {
            write!(abaqus_input_buffer, "{},{},{},{}\n",i+1,node_vec[i][0],node_vec[i][1],node_vec[i][2]);
        }
        write!(abaqus_input_buffer,"*Element, type=DC3D8\n");
        for i in 0..elem_vec.len(){
            write!(abaqus_input_buffer,"{},{},{},{},{},{},{},{},{}\n",i+1, elem_vec[i][0]+1,elem_vec[i][1]+1,elem_vec[i][3]+1,elem_vec[i][2]+1,elem_vec[i][4]+1,elem_vec[i][5]+1,elem_vec[i][7]+1,elem_vec[i][6]+1);
        }

        write!(abaqus_input_buffer,
               "*Nset, nset=Set-1, generate
    1,  {},     1
 *Elset, elset=Set-1, generate
  1,  {},   1",node_vec.len(),elem_vec.len());
        write!(abaqus_input_buffer,
               "
** Section: Section-1
*Solid Section, elset=Set-1, material=Material-1
,
*End Instance
**
*ELEMENT PROGRESSIVE ACTIVATION, NAME = \"EPA\", ELSET = PART-1-1.SET-1
*End Assembly
**
** MATERIALS
**
*Material, name=Material-1
*Density
1240,
*Conductivity
 0.33,
*Specific Heat
2020.,
");


        write!(abaqus_input_buffer,"
** ----------------------------------------------------------------
**
** STEP: Step-1
**
*Step, name=Step-1, nlgeom=NO, inc=1000
*Heat Transfer, end=PERIOD
1., 100, , ,
",
        );


    }

    pub fn write_abaqus_input_am_file(&self, abaqus_meshfile: &str, activated_elements:Vec<usize>, activation_times_only: Vec<f64>, cooldown_period:f64, timestep:f64){
        let mut abaqus_input = File::create(abaqus_meshfile).expect("can't create the nodefile at the specified location");
        let mut abaqus_input_buffer = BufWriter::with_capacity(100000,abaqus_input);

        let mut node_table_old_to_new = HashMap::with_capacity(100000);
        let mut node_vec = Vec::with_capacity(100000);
        let mut elem_vec = Vec::with_capacity(100000);
        let mut node_counter = 0;

        for i in activated_elements{
            let elem = self.celllist[i];
            let mut elem_new:[usize;8] = [0;8];
            for j in 0..8{
                if node_table_old_to_new.contains_key(&elem[j]){
                    let new_node_number = node_table_old_to_new.get(&elem[j]).expect("cant find key");
                    elem_new[j] = *new_node_number;
                } else{
                    node_table_old_to_new.insert(elem[j].clone(),node_counter.clone());
                    elem_new[j] = node_counter.clone();
                    node_vec.push([self.nodelist[elem[j]].pt.x,self.nodelist[elem[j]].pt.y,self.nodelist[elem[j]].pt.z]);
                    node_counter += 1;
                }
            }
            elem_vec.push(elem_new);
        }

        write!(abaqus_input_buffer,
               "*Heading
** Job name: abaqus_meshfile Model name: Model-1
** Generated by: Sunil Bhandari
*Preprint, echo=NO, model=NO, history=NO, contact=NO
**
** PARTS
**
*Part, name=Part-1
*End Part
**
**
** ASSEMBLY
**
*Assembly, name=Assembly
**
*Instance, name=Part-1-1, part=Part-1
*Node\n");
        for i in 0..node_vec.len() {
            write!(abaqus_input_buffer, "{},{},{},{}\n",i+1,node_vec[i][0],node_vec[i][1],node_vec[i][2]);
        }
        write!(abaqus_input_buffer,"*Element, type=DC3D8\n");
        for i in 0..elem_vec.len(){
            write!(abaqus_input_buffer,"{},{},{},{},{},{},{},{},{}\n",i+1, elem_vec[i][0]+1,elem_vec[i][1]+1,elem_vec[i][3]+1,elem_vec[i][2]+1,elem_vec[i][4]+1,elem_vec[i][5]+1,elem_vec[i][7]+1,elem_vec[i][6]+1);
        }

        write!(abaqus_input_buffer,
               "*Nset, nset=Set-1, generate
    1,  {},     1
*Elset, elset=Set-1, generate
  1,  {},   1
** Section: Section-1
*Solid Section, elset=Set-1, material=Material-1
,
*End Instance
**
",node_vec.len(),elem_vec.len());
        write!(abaqus_input_buffer,
               "
*ELEMENT PROGRESSIVE ACTIVATION, NAME = \"EPA\", ELSET = PART-1-1.SET-1, FOLLOW DEFORMATION=YES
*End Assembly
**
** MATERIALS
**
*Material, name=Material-1
*Density
1240,
*Conductivity
 0.33,
*Specific Heat
2020.,
**PROPERTY TABLE DEFINITIONS
*PARAMETER TABLE TYPE, NAME = \"ACTIVATION SERIES\", PARAMETERS = 1
FLOAT,,,\"Activation Time\",
*TABLE COLLECTION, NAME = \"ACTIVATION TABLE\"
*PARAMETER TABLE, TYPE = \"ACTIVATION SERIES\" \n");

        for i in 0..activation_times_only.len(){
            write!(abaqus_input_buffer, "{},\n", activation_times_only[i]);
        }
        write!(abaqus_input_buffer,"
** PREDEFINED FIELDS
**
** Name: Predefined Field-1   Type: Temperature
*Initial Conditions, type=TEMPERATURE
Part-1-1.Set-1, 220.
** ----------------------------------------------------------------
**
** STEP: Step-1
**
*Step, name=Step-1, nlgeom=NO, inc={}
*Heat Transfer, end=PERIOD
{}, {}, , ,
*ACTIVATE ELEMENTS, ACTIVATION = \"EPA\"
\"ACTIVATION TABLE\"
*FILM
Part-1-1.Set-1, FFS,25,100
**
**
** OUTPUT REQUESTS
**
*Restart, write, frequency=0
**
** FIELD OUTPUT: F-Output-1
**
*Output, field
*Node Output
NT, RFL
*Element Output, directions=YES
EACTIVE, HFL, TEMP
**
** HISTORY OUTPUT: H-Output-1
**
*Output, history
*Radiation Output
FTEMP,
*End Step",
               activation_times_only.len()+activation_times_only.len()/10, timestep,
               activation_times_only[activation_times_only.len()-1]+cooldown_period);


    }

    pub fn write_abaqus_UEPACTIVATIONVOL_file(&self, uepactivationvolfile: &str, num_act_elems: usize){
        let mut uepactvol = File::create(uepactivationvolfile).expect("can't create the nodefile at the specified location");
        let mut uepactvolbuf = BufWriter::with_capacity(100000,uepactvol);
        write!(uepactvolbuf,
               "      subroutine uepactivationvol(
     * lFlags,
     * epaName,
     * noel,
     * nElemNodes,
     * iElemNodes,
     * mcrd,
     * coordNodes,
     * uNodes,
     * kstep,
     * kinc,
     * time,
     * dtime,
     * temp,
     * npredef,
     * predef,
     * nsvars,
     * svars,
     * sol,
     * solinc,
     * volFract,
     * nVolumeAddEvents,
     * volFractAdded,
     * csiAdded)
      include 'aba_param.inc'
	  include 'aba_tcs_param.inc'
      dimension
     * lFlags(*),
     * iElemNodes(nElemNodes),
     * coordNodes(mcrd,nElemNodes),
     * uNodes(mcrd,nElemNodes),
     * time(2),
     * temp(2,nElemNodes),
     * predef(2,npredef,nElemNodes),
     * svars(2,nsvars),
     * sol(nElemNodes),
     * solinc(nElemNodes),
     * volFract(*),
     * volFractAdded(*),
     * csiAdded(3,*)
	  character*80 parameterTableLabel, cParams(1)
      dimension iParamsDataType(1), iParams(1), rParams({})
	  t1 = time(1)
	  if (noel.lt.{}) then
	  call getParameterTableRow(\'ACTIVATION SERIES\', noel, numParams, iParamsDataType, iParams,
     * rParams, cParams, jError)
	  if ((rParams(1).lt.t1).AND.(volFract(1).lt.1.0)) then
		volFractAdded(1) = 1.0;
	  end if
	  end if
      return
      end", num_act_elems, num_act_elems  );

    }


}