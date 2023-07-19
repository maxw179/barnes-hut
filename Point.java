import java.util.ArrayList;
import java.io.Serializable;

public class Point
{
  private double[] x;
  private double[] v;
  private double m;

  public Point(double[] x, double[] v, double m)
  {
    this.x = x;
    this.v = v;
    this.m = m;
  }

  public Point(Point p)
  {
    this.x = new double[2];
    this.v = new double[2];
    this.m = p.m;
    for(int i = 0; i < 2; i++)
    {
      x[i] = p.x[i];
      v[i] = p.v[i];
    }
  }

  public void update(double[] f, double dt)
  {
    for(int i = 0; i < 2; i++)
    {
      v[i] = v[i] + (f[i]/m) * dt;
      x[i] = x[i] + v[i] * dt + .5 * (f[i]/m) * dt * dt;
    }
  }

  public double[] getX()
  {
    return x;
  }
  public double[] getV()
  {
    return v;
  }
  public double getM()
  {
    return m;
  }

  public String toString()
  {
    return "[" + x[0] + ", " + x[1] + "]";
  }

}
