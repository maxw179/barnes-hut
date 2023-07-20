import java.util.ArrayList;
import java.lang.Math;
import java.io.*;

public class Runner
{
  public static void main(String[] args)
  {
    if(args.length != 0)
    {
      //Initialize Simulation Parameters
      int nBodies = 1000;
      double xScalar = 0.7;
      double mass = 5;
      int timeSteps = 1000;
      double dt = 1.5;
      double theta = 1.5;
      //Try to save the simulation to _filename.csv
      runSimulation(nBodies, xScalar, mass, timeSteps, dt, theta, "_" + args[0] + ".csv");

    }
    else
    {
      System.out.println("Please include filename.");
    }
  }

  public static void runSimulation(int nBodies, double xScalar, double mass, int timeSteps, double dt, double theta, String fileName)
  {
    try
    {
      //initialize FileWriter
      FileWriter writer = new FileWriter(fileName);
      //initialize points
      //ArrayList<Point> points = generateGalaxy(nBodies, 0.5, 0.5, xScalar, mass);
      ArrayList<Point> points = generateRandomSquare(nBodies, xScalar, mass);
      long startTime = System.nanoTime();
      long endTime = System.nanoTime();
      double elapsedSeconds = (endTime - startTime) * 1e-9;
      System.out.printf("Initialization completed in: %.3fs%n", elapsedSeconds);
      //write initial points to the file
      writeCSV(points, writer);

      //loop through steps of the simulation
      for(int i = 0; i < timeSteps; i++)
      {
        //run a step of the simulation
        startTime = System.nanoTime();
        ArrayList<Point> newPoints = oneStep(points, theta, dt);
        endTime = System.nanoTime();
        elapsedSeconds = (endTime - startTime) * 1e-9;
        System.out.printf("Step " + (i + 1) + " completed in: %.3fs%n", elapsedSeconds);

        //add the new points to the file
        writeCSV(newPoints, writer);
        //set the old points to be the new points
        points = newPoints;
        System.out.println(points.size() + " Bodies remain.");
      }
      writer.close();
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }

  }

  public static void writeCSV(ArrayList<Point> points, FileWriter writer)
  {
    try
    {
      for(Point point : points)
      {
        String x = String.valueOf(point.getX()[0]).substring(0, 5) + ", " + String.valueOf(point.getX()[1]).substring(0, 5) + ", ";
        writer.append(x);
      }
      writer.append('\n');
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
  }

  public static ArrayList<Point> oneStep(ArrayList<Point> points, double theta, double dt)
  {
    //find out how large the quadtree should be
    double[] minsMaxes = getMinsMaxes(points);
    double xMin = minsMaxes[0];
    double xMax = minsMaxes[1];
    double yMin = minsMaxes[2];
    double yMax = minsMaxes[3];

    //create a new quadtree and add the old points into it
    QuadTree q = new QuadTree(max(xMax - xMin, yMax - yMin));
    q.add(points, xMin, yMin);
    //get the points and forces from the quadtree
    ArrayList<Point> reformattedPoints = q.getPoints();
    ArrayList<double[]> forces = q.getForces(theta);

    //new arraylist to store the updated points
    ArrayList<Point> newPoints = new ArrayList<Point>();
    //iterate through each point
    for(int i = 0; i < reformattedPoints.size(); i++)
    {
      //make a deep copy of the old point
      Point oldPoint = reformattedPoints.get(i);
      Point newPoint = new Point(oldPoint);
      //update the new point
      newPoint.update(forces.get(i), dt);
      //add the new point into new points
      newPoints.add(newPoint);
    }
    return newPoints;
  }

  //return the max of two doubles
  public static double max(double a, double b)
  {
    if(a > b)
    {
      return a;
    }
    return b;
  }

  //returns the xmin, xmax, ymin, and ymax of a list of points
  public static double[] getMinsMaxes(ArrayList<Point> points)
  {
    double xMin = points.get(0).getX()[0];
    double xMax = points.get(0).getX()[0];
    double yMin = points.get(0).getX()[1];
    double yMax = points.get(0).getX()[1];
    //iterate over all the points to identify the x min, x max, y min, and y max
    for(Point point : points)
    {
      if(point.getX()[0] < xMin)
      {
        xMin = point.getX()[0];
      }
      if(point.getX()[0] > xMax)
      {
        xMax = point.getX()[0];
      }
      if(point.getX()[1] < yMin)
      {
        yMin = point.getX()[1];
      }
      if(point.getX()[1] > yMax)
      {
        yMax = point.getX()[1];
      }
    }
    double[] minsMaxes = {xMin, xMax, yMin, yMax};
    return minsMaxes;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  public static ArrayList<Point> generateRandomSquare(int nBodies, double xScalar, double mass)
  {
    //create a new list of points
    ArrayList<Point> list = new ArrayList<Point>();
    //iterate through each point
    for(int i = 0; i < nBodies; i++)
    {
      //initialize the point with an x, y, and velocity
      double[] x = generateX(xScalar);
      Point p = new Point(x, generateV(0.0), mass);
      list.add(p);
    }
    //return the list
    return list;
  }

  public static ArrayList<Point> generateBigBang(int nBodies, double xScalar, double mass)
  {
    //create a new list of points
    ArrayList<Point> list = new ArrayList<Point>();
    double pi = 3.14159265359;
    //iterate through each point
    for(int i = 0; i < nBodies; i++)
    {
      //initialize the point with an x, y, and velocity
      double theta = Math.random() * (2 * pi);
      double rRandom = Math.random() * 0.01 - 0.005;
      double[] x = {(xScalar + rRandom) * Math.cos(theta) + .5, (xScalar + rRandom) * Math.sin(theta) + .5};
      double xV = -(x[0] - .5);
      double yV = -(x[1] - .5);
      double[] v = {xV/10.0, yV/10.0};
      Point p = new Point(x, v, mass);
      list.add(p);
    }
    //return the list
    return list;
  }

  public static double[] circularX(double xScalar)
  {
    double[] pair = {1.0,1.0};
    double offset = (.5);
    while(Math.sqrt(pair[0] * pair[0] + pair[1] * pair[1]) > xScalar)
    {
      pair[0] = Math.random() * (2 * xScalar) - xScalar;
      pair[1] = Math.random() * (2 * xScalar) - xScalar;
    }
    pair[0] += offset;
    pair[1] += offset;
    return pair;
  }

  public static double[] generateX(double xScalar)
  {
    double x = Math.random() * xScalar + (.5 - xScalar/2);
    double y = Math.random() * xScalar + (.5 - xScalar/2);
    double[] pair = {x, y};
    return pair;
  }

  public static double[] generateV(double xScalar)
  {
    double x = Math.random() * (xScalar) + (-xScalar/2);
    double y = Math.random() * (xScalar) + (-xScalar/2);
    double[] pair = {x, y};
    return pair;
  }

  public static ArrayList<Point> generateSolarSystem(int nPlanets, double xScalar, double solarMass)
  {
    //create a new list of points
    ArrayList<Point> list = new ArrayList<Point>();
    //generate the sun at the center of the simulation
    double[] sunX = {xScalar/2, xScalar/2};
    double[] sunV = {0.0,0.0};
    list.add(new Point(sunX, sunV, solarMass));

    //generate planets in a line
    for(int i = 0; i < nPlanets; i++)
    {
      //calculate speed according to sqrt(GM/r)
      double[] px = {xScalar/2 + xScalar/8 + (xScalar/3) * ((double)i/nPlanets), xScalar/2};
      double[] pv = {0.0, Math.sqrt((6.67e-11 * solarMass)/(xScalar/8 + xScalar/3 * ((double)i/nPlanets)))};
      //infinitessimal mass so no cross-planet gravitation
      double planetMass = 1.0e-25;
      list.add(new Point(px, pv, planetMass));
    }
    return list;
  }

  public static ArrayList<Point> generateGalaxy(int nBodies, double xOffset, double yOffset, double xScalar, double mass)
  {
    ArrayList<Point> list = new ArrayList<Point>();
    double pi = 3.14159265359;
    for(int i = 0; i < nBodies; i++)
    {
      //random r and theta
      double theta = Math.random() * (2 * pi);
      double lambda = 4.0;
      double r = getExponentialRandom(lambda) * (xScalar/2);
      //get x and y from theta
      double x = r * Math.cos(theta) + (xOffset);
      double y = r * Math.sin(theta) + (yOffset);
      double[] X = {x, y};

      //find mass contained within the radius
      //(integral of probability from 0 to r) divided by (integral of probability from 0 to xScalar/2) all times mass
      double massContained = (integralExponential(lambda, r) / integralExponential(lambda, xScalar/2)) * (mass * nBodies);
      double scalarVelocity = Math.sqrt((6.67e-11 * massContained)/r);
      //calculate velocity according to contained mass
      double[] V = {Math.sin(theta) * scalarVelocity, -Math.cos(theta) * scalarVelocity};
      list.add(new Point(X, V, mass));
    }
    return list;
  }

  public static double integralExponential(double lambda, double x)
  {
    return (-1/lambda) * Math.exp(-lambda * x) + (1/lambda);
  }

  public static double getExponentialRandom(double lambda)
  {
    return Math.exp(-lambda * Math.random());
  }


}
