import java.util.ArrayList;
import java.lang.Math;
import java.io.*;

public class Node
{
  private Point point;
  private double width;
  private Node[] children;
  private double[] centerOfMass;
  private double totalMass;

  public Node(double width, Point point)
  {
    this.width = width;
    this.point = point;
    this.centerOfMass = new double[2];
    this.totalMass = 0.0;
    children = new Node[4];
  }

  public Node(double width)
  {
    this(width, null);
  }

  //CLASS METHODS
  public boolean isExternal()
  {
    return point != null;
  }

  public double getWidth()
  {
    return width;
  }

  public Point getPoint()
  {
    return point;
  }

  public void setChild(Node newChild, int index)
  {
    children[index] = newChild;
  }

  public Node[] getChildren()
  {
    return children;
  }

  public Node getChild(int index)
  {
    return children[index];
  }

  public double getMass()
  {
    //special double checking
    if(totalMass > 0.00001)
    {
      return totalMass;
    }
    calculateMass();
    return totalMass;
  }

  public double[] getCOM()
  {
    //special double checking
    if(centerOfMass[0] > 0.00001)
    {
      return centerOfMass;
    }
    calculateCOM();
    return centerOfMass;
  }

  public ArrayList<Point> getPoints()
  {
    ArrayList<Point> points = new ArrayList<Point>();
    if(isExternal())
    {
      points.add(point);
    }
    else
    {
      for(Node child : children)
      {
        if(child != null)
        {
          points.addAll(child.getPoints());
        }
      }
    }
    return points;
  }

  public double[] getDistanceVector(Node other)
  {
    double[] thisPos = getCOM();
    double[] otherPos = other.getCOM();

    double[] distance = {otherPos[0] - thisPos[0], otherPos[1] - thisPos[1]};
    return distance;
  }

  public double getDistanceScalar(double[] distanceVector)
  {
    return Math.sqrt(Math.pow(distanceVector[0], 2) + Math.pow(distanceVector[1], 2));
  }

  public double getDistanceScalar(Node other)
  {
    return getDistanceScalar(getDistanceVector(other));
  }

  public double[] getForce(Node other)
  {
    double[] distanceVector = getDistanceVector(other);
    double distanceScalar = getDistanceScalar(distanceVector);
    //avoid divison by 0 distance
    if(distanceScalar < 0.000001)
    {
      return new double[2];
    }
    //normal case
    else
    {
      double massProduct = getMass() * other.getMass();
      double forceScalar = (6.67e-11 * massProduct)/Math.pow(Math.pow(distanceScalar, 2) + .001, 1.5);
      double[] forceVector = {forceScalar * distanceVector[0], forceScalar * distanceVector[1]};
      return forceVector;
    }
  }

  public double[] recursiveForce(Node currentNode, double theta)
  {
    double sdRatio = currentNode.getWidth() / getDistanceScalar(currentNode);
    if(currentNode.isExternal() || sdRatio < theta)
    {
      return getForce(currentNode);
    }
    else
    {
      double xForce = 0.0;
      double yForce = 0.0;
      for(Node child : currentNode.getChildren())
      {
        if(child != null)
        {
          double[] childForce = recursiveForce(child, theta);
          xForce += childForce[0];
          yForce += childForce[1];
        }
      }
      double[] totalForce = {xForce, yForce};
      return totalForce;
    }
  }

  public ArrayList<double[]> getForces(Node root, double theta)
  {
    ArrayList<double[]> forces = new ArrayList<double[]>();
    if(isExternal())
    {
      forces.add(recursiveForce(root, theta));
    }
    else
    {
      for(Node child : children)
      {
        if(child != null)
        {
          forces.addAll(child.getForces(root, theta));
        }

      }
    }
    return forces;
  }

  public void calculateMass()
  {
    if(isExternal())
    {
      totalMass = point.getM();
    }
    else
    {
      double sum = 0.0;
      for(Node child : children)
      {
        if(child != null)
        {
          sum += child.getMass();
        }
      }
      totalMass = sum;
    }
  }

  public void calculateCOM()
  {
    if(isExternal())
    {
      centerOfMass = point.getX();
    }
    else
    {
      double x = 0.0;
      double y = 0.0;
      for(Node child : children)
      {
        if(child != null)
        {
          double mass = child.getMass();
          x += child.getCOM()[0] * mass;
          y += child.getCOM()[1] * mass;
        }
      }
      double[] result = {x/getMass(), y/getMass()};
      centerOfMass = result;
    }
  }

  public String toString()
  {
    if(isExternal())
    {
      return "[" + point.toString() + "]";
    }
    else
    {
      String s = "(";
      for(Node child : children)
      {
        if(child != null)
        {
          s += child.toString();
        }
        s += ", ";
      }
      return s + ")";
    }
  }
}
