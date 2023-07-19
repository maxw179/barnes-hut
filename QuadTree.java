import java.util.ArrayList;

public class QuadTree
{
  private Node root;

  public QuadTree(double width)
  {
    root = new Node(width);
  }

  public String toString()
  {
    return root.toString();
  }

  public void add(ArrayList<Point> points)
  {
    for(Point point : points)
    {
      add(point, root, 0.0, 0.0);
    }
  }

  public ArrayList<Point> getPoints()
  {
    return root.getPoints();
  }

  public ArrayList<double[]> getForces(double theta)
  {
    return root.getForces(root, theta);
  }

  private void add(Point point, Node last, double xOffset, double yOffset)
  {
    double xValue = point.getX()[0];
    double yValue = point.getX()[1];
    double width = last.getWidth();

    boolean[] xConditionTable = {(xValue <= xOffset + width/2), (xValue > xOffset + width/2), (xValue <= xOffset + width/2), (xValue > xOffset + width/2)};
    boolean[] yConditionTable = {(yValue > yOffset + width/2), (yValue > yOffset + width/2), (yValue <= yOffset + width/2), (yValue <= yOffset + width/2)};

    double[] xOffsetTable = {xOffset, xOffset + width/2, xOffset, xOffset+ width/2};
    double[] yOffsetTable = {yOffset + width/2, yOffset + width/2, yOffset, yOffset};

    String[] printTable = {"NORTHWEST", "NORTHEAST", "SOUTHWEST", "SOUTHEAST"};

    for(int i = 0; i < 4; i++)
    {
      if(xConditionTable[i] && yConditionTable[i])
      {
        //System.out.println(printTable[i]);
        //case if open spot
        if(last.getChild(i) == null)
        {
          last.setChild((new Node(width/2, point)), i);
        }
        //case if occupied spot
        else
        {
          //case if the occupied space is external
          if(last.getChild(i).isExternal())
          {
            Point collisionPoint = last.getChild(i).getPoint();
            last.setChild((new Node(width/2)), i);
            add(collisionPoint, last.getChild(i), xOffsetTable[i], yOffsetTable[i]);
          }
          //from here, treat the occupied space as internal
          add(point, last.getChild(i), xOffsetTable[i], yOffsetTable[i]);
        }

      }
    }
  }
}
