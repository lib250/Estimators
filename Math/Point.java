public class Point{

    private double x;
    private double y;
    private double z;

    public Point(double x, double y){
        this.x = x;
        this.y = y;
        this.z = 0;
    }

    public Point(double x, double y, double z){
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public double x(){
        return x;
    }

    public double y(){
        return y;
    }

    public double z(){
        return z;
    }

    public double distanceTo(Point other){
        return Math.sqrt((other.x() - x) * (other.x() - x) +
                         (other.y() - y) * (other.y() - y) +
                         (other.z() - z) * (other.z() - z));
    }
}