use crate::primitives::Point;
use std::fs::File;
use std::io::{BufReader, BufRead, BufWriter};
use std::iter::FromIterator;
use std::io::Write;

/// GCodeReader stores the information from GCode file in form of move segments, is extrusion on for
/// the segments, and speed of movement for those segments.
#[derive(Debug,Clone)]
pub(crate) struct GCodeReader{
    pub(crate) segment: Vec<[Point;2]>,
    pub(crate) speed: Vec<f64>,
    pub(crate) is_extrusion_on: Vec<bool>,
    pub(crate) xmin:f64,
    pub(crate) xmax:f64,
    pub(crate) ymin:f64,
    pub(crate) ymax:f64,
    pub(crate) zmin:f64,
    pub(crate) zmax:f64,
}

/// implements methods for GCodeReader
impl GCodeReader{
    /// create new GCodeReader by reading information from a GCodeFile. Relevant information needs
    /// to start with line ";Printing starts here" and end with line "; layer end"
    pub fn new(filename: &str)-> GCodeReader{
        let file_output = GCodeReader::read_file(filename);

        //let zz = 0.0;
        let mut segment = Vec::new();
        let mut speed = Vec::new();
        let mut is_extrusion_on = Vec::new();
        let mut Point1 = Point{x:0.0,y:0.0,z:0.0,};
        let mut Point2 = Point{x:0.0,y:0.0,z:0.0,};
        Point1 = file_output.0[1].clone();
        let mut xmin = f64::INFINITY;  let mut xmax=f64::NEG_INFINITY;
        let mut ymin= f64::INFINITY;  let mut ymax =f64::NEG_INFINITY;
        let mut zmin= f64::INFINITY;  let mut zmax =f64::NEG_INFINITY;
        for i in 2..file_output.0.len(){
            Point2 = file_output.0[i].clone();
            let is_true = file_output.2[i];
            if file_output.2[i] {
                if Point2.x > xmax { xmax = Point2.x };
                if Point2.y > ymax { ymax = Point2.y };
                if Point2.z > zmax { zmax = Point2.z }
                if Point2.x < xmin { xmin = Point2.x };
                if Point2.y < ymin { ymin = Point2.y };
                if Point2.z < zmin { zmin = Point2.z }
            }
            segment.push([Point1.clone(),Point2.clone()]);
            Point1 = Point2;
            speed.push(file_output.1[i]);
            is_extrusion_on.push(file_output.2[i]);
        }

        return GCodeReader{
            segment,
            speed,
            is_extrusion_on,
            xmin,
            xmax,
            ymin,
            ymax,
            zmin,
            zmax,
        }
    }

    /// reads the GCode file and gets the movement path, movement speed, extrusion on/off, extrusion
    /// distance
    pub fn read_file(filename:&str)->(Vec<Point>,Vec<f64>,Vec<bool>,Vec<f64>){
        println!("gcode filename -> {}", filename);
        let file = File::open(filename).expect("cant open the file");
        let bfile = BufReader::with_capacity(10000,&file);
        let mut move_vector = Vec::with_capacity(10000);
        let mut speed_vector = Vec::with_capacity(10000);
        let mut is_extrusion_on_vector = Vec::with_capacity(10000);
        let mut extrusion_vector = Vec::with_capacity(10000);
        let mut start = false;
        let mut extrude = false;
        let mut movespeed:f64 = 0.0;
        let mut x1:f64 = 0.0;
        let mut y1:f64 = 0.0;
        let mut z1:f64 = 0.0;
        let mut extrusion_length :f64 = 0.0;
        let mut moveupdate = false;

        for i in bfile.lines(){
            let k = i.unwrap();
            if k== ";Printing starts here".to_string() {start = true;}
            if k== "; layer end".to_string() {start = false;}
            extrude = false;
            moveupdate = false;


            if start {
                let mut m = k.split_whitespace();

                let kk = m.next().unwrap();
                if kk == "G1".to_string() {
                    moveupdate = false;
                    let mut mclone = m.clone();
                    let k2clone = mclone.next().unwrap();
                    let mut k2chars = k2clone.chars();
                    let nextk2chars = k2chars.next().unwrap();
                    if (nextk2chars == 'X' || nextk2chars == 'Y' || nextk2chars =='Z') {
                        moveupdate = true;

                        for k2 in m {
                            // let k2 = m.next().unwrap();
                            let mut kchars = k2.chars();
                            let nextchar = kchars.next().unwrap();

                            if nextchar == 'Z' {
                                z1 = String::from_iter(kchars).parse::<f64>().unwrap();
                            } else if nextchar == 'X' {
                                x1 = String::from_iter(kchars).parse::<f64>().unwrap();
                            } else if nextchar == 'Y' {
                                match String::from_iter(kchars).parse::<f64>() {
                                    Err(E) => {
                                        println!("err.");
                                        //let zzz = 000;
                                    }

                                    Ok(value) => y1 = value,
                                }
                            } else if nextchar == 'E' {
                                extrude = true;
                                extrusion_length = String::from_iter(kchars).parse::<f64>().unwrap();
                            } else if nextchar == 'F' {
                                movespeed = String::from_iter(kchars).parse::<f64>().unwrap();
                            }
                        }
                    }
                    if moveupdate {
                        move_vector.push(Point { x: x1.clone()*1e-3, y: y1.clone()*1e-3, z: z1.clone()*1e-3 });
                        speed_vector.push(movespeed.clone()/60.0*1e-3);
                        is_extrusion_on_vector.push(extrude);
                        extrusion_vector.push(extrusion_length*1e-3);
                    }
                }
            }

        }
        return (move_vector,speed_vector,is_extrusion_on_vector,extrusion_vector);

    }

    pub fn write_abaqus_input(&self,filename:&str){
        let mut codedump = File::create(filename).expect("can't create the specified file");
        let mut codedumpbuffer = BufWriter::with_capacity(100000,codedump);
        let mut totaltime:f64 = 0.0;
        for i in 0..self.segment.len(){
            let time = self.segment[i][0].dist_to(&self.segment[i][1])/self.speed[i];

            if (self.is_extrusion_on[i]) {
                write!(codedumpbuffer,"{},{},{},{},{}\n",totaltime, self.segment[i][0].x+0.0004,
                       self.segment[i][0].y,self.segment[i][0].z,1);
            }
            else{
                write!(codedumpbuffer,"{},{},{},{},{}\n",totaltime, self.segment[i][0].x+0.0004,
                       self.segment[i][0].y,self.segment[i][0].z,0);
            }

            totaltime = totaltime + time;
        }
        write!(codedumpbuffer,"{},{},{},{},{}\n",totaltime, self.segment[self.segment.len()-1][1].x+0.0004,
               self.segment[self.segment.len()-1][1].y,
               self.segment[self.segment.len()-1][1].z,1);
    }
}