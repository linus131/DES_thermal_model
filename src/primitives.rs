use std::fmt;

/// Point struct holds x, y, and z coordinates as f64
#[derive(Debug, Copy, Clone)]
pub struct Point{
    pub x:f64,
    pub y:f64,
    pub z:f64,
}
/// Methods for Point struct
impl Point{
    ///distance from one point to another
    pub fn dist_to(&self, other:&Point)->f64{
        let xd = self.x - other.x;
        let yd = self.y - other.y;
        let zd = self.z - other.z;
        return (xd*xd+yd*yd+zd*zd).sqrt();
    }
}


impl fmt::Display for Point{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {}, {})", self.x, self.y, self.z)
    }
}


/// Node has a point and and an index
#[derive(Debug, Copy, Clone)]
pub struct Node{
    pub(crate) pt:Point,
    pub(crate) index:usize,
}

/// Cell has unique index, eight nodes,
/// material properties, and the center point of the cell
#[derive(Debug,Copy,Clone)]
struct Cell{
    index:usize,
    sp_heat_cap:f64,
    volume:f64,
    density:f64,
    kx:f64,
    kyz:f64,
    center:Point,
    nodes:[usize;8],
}
