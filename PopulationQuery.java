import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Scanner;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.*;

//VERSION 2 start
//this class wraps the infomation of the four corners so that we can pass it
//during the recursion
class FourCorner{
	float top;
	float bottom;
	float left;
	float right;
	public FourCorner(float top, float bottom, float left, float right){
		this.top = top;
		this.bottom = bottom;
		this.left = left;
		this.right = right;
	}
}
//we use this class to use the Fork-join framework to parallelize the computation
//of the four corners
class V2FindCorner extends RecursiveTask<FourCorner>{
	private static final int CUTOFF = 100;
	private int low;
	private int high;
	private CensusGroup[] arr;
	private int dataSize;
	public V2FindCorner(int low, int high, CensusGroup[] arr, int dataSize){
		this.low = low;
		this.high = high;
		this.arr = arr;
		this.dataSize = dataSize;
	}
	@Override
	public FourCorner compute(){
		//base case
		if (high - low <= CUTOFF - 1){
			float top = Integer.MIN_VALUE;
			float bottom = Integer.MAX_VALUE;
			float left = Integer.MAX_VALUE;
			float right = Integer.MIN_VALUE;

			for (int i = low; i <= high; i++){
				CensusGroup group = this.arr[i];
				//System.out.println(group.latitude);
				top = Math.max(top, group.latitude);
				bottom = Math.min(bottom, group.latitude);
				left = Math.min(left, group.longitude);
				right = Math.max(right, group.longitude);
			}
			FourCorner fourCorner = new FourCorner(top, bottom, left, right);
			return fourCorner;
		}
		//we recursively calc the four corners of the US
		int mid = (low + high)/2;
		//System.out.println(mid);
		V2FindCorner leftThread = new V2FindCorner(low, mid, arr, dataSize);
		V2FindCorner rightThread = new V2FindCorner(mid + 1, high, arr, dataSize);

		leftThread.fork();
		FourCorner rightFourCorner = rightThread.compute();
		FourCorner leftFourCorner = leftThread.join();

		float top = Math.max(rightFourCorner.top, leftFourCorner.top);
		float bottom = Math.min(rightFourCorner.bottom, leftFourCorner.bottom);
		float left = Math.min(rightFourCorner.left, leftFourCorner.left);
		float right = Math.max(rightFourCorner.right, rightFourCorner.left);
		return new FourCorner(top, bottom, left, right);
	}
	//return the calculated corners
	public float[] findCorner2() throws InterruptedException{
		ForkJoinPool fjPool = new ForkJoinPool();
		FourCorner result = fjPool.invoke(new V2FindCorner(0, dataSize - 1, arr, dataSize));
		float[] corners = new float[4];
		corners[0] = result.top;
		corners[1] = result.bottom;
		corners[2] = result.left;
		corners[3] = result.right;
		return corners;
	}
}
//this class calculate the distribution of the population based on blocks we assigned
//and the total population
class V2FindPopulation extends RecursiveTask<Pair<Integer, int[][]>>{
	private static final int CUTOFF = 100;
	private int low;
	private int high;
	private CensusGroup[] arr;
	private int dataSize;
	private float[] corners;
	private int rows;
	private int cols;

	public V2FindPopulation(int low, int high, CensusGroup[] arr, int dataSize, float[] corners, int cols, int rows){
		this.low = low;
		this.high = high;
		this.arr = arr;
		this.dataSize = dataSize;
		//this.matrixPeople =
		this.corners = corners;
		this.rows = rows;
		this.cols = cols;
	}
	@Override
	public Pair<Integer, int[][]> compute(){
		//base case
		if (high - low <= CUTOFF - 1){
			int totalPopulation = 0;
			int[][] matrixPeople = new int[rows][cols];
			for(int i = low; i <= high; i++) {
				CensusGroup group = arr[i];
				int[] index = findIndex(group.longitude, group.latitude, corners, cols, rows);
				int rowNum = index[0];
				int colNum = index[1];
				if (rowNum == rows && colNum == cols){
					matrixPeople[rowNum - 1][colNum - 1] += group.population;
				}else if (rowNum == rows){
					matrixPeople[rowNum - 1][colNum] += group.population;
				}else if (colNum == cols){
					matrixPeople[rowNum][colNum - 1] += group.population;
				}else{
					matrixPeople[rowNum][colNum] += group.population;
				}
				totalPopulation += group.population;
			}
			return new Pair<Integer, int[][]>(totalPopulation, matrixPeople);
		}
		//divide
		int mid = (low + high)/2;
		V2FindPopulation leftThread = new V2FindPopulation(low, mid, arr, dataSize, corners, cols, rows);
		V2FindPopulation rightThread = new V2FindPopulation(mid + 1, high, arr, dataSize, corners, cols, rows);
		//conquer
		leftThread.fork();
		Pair<Integer, int[][]> rightPair = rightThread.compute();
		Pair<Integer, int[][]> leftPair = leftThread.join();
		int result = rightPair.getElementA() + leftPair.getElementA();
		int[][] matrixPeople = leftPair.getElementB();
		for (int i = 0; i < matrixPeople.length; i++){
			for (int j = 0; j < matrixPeople[0].length; j++){
				matrixPeople[i][j] += rightPair.getElementB()[i][j];
			}
		}

		return new Pair<Integer, int[][]>(result, matrixPeople);
	}
	//return the calculated total population and the distribution of the population
	public Pair<Integer, int[][]> findPopulation2() throws InterruptedException{
		ForkJoinPool fjPool = new ForkJoinPool();
		return fjPool.invoke(new V2FindPopulation(0, dataSize - 1, arr, dataSize, corners, cols, rows));
	}
	//find the index in matrix
	public int[] findIndex(float longitude, float latitude, float[] corners, int cols, int rows) {
		float top = corners[0];
		float bottom = corners[1];
		float left = corners[2];
		float right = corners[3];

		float height = (top - bottom)/rows;
		float width = (right - left)/cols;

		int[] index = new int[2];

		index[0] = (int)((latitude-bottom)/height);


		index[1] = (int)((longitude-left)/width);


		return index;
	}
}

class V5FindPopulation extends Thread{
	private int low;
	private int high;
	private CensusGroup[] arr;
	private int dataSize;
	private float[] corners;
	private int rows;
	private int cols;
	public static int[][] matrixPeople;
	private static int totalPopulation = 0;
  public static Object lock = new Object();

	public V5FindPopulation(int low, int high, CensusGroup[] arr, int dataSize, float[] corners, int cols, int rows){
		this.low = low;
		this.high = high;
		this.arr = arr;
		this.dataSize = dataSize;
		//this.matrixPeople =
		this.corners = corners;
		this.rows = rows;
		this.cols = cols;

	}
	public Pair<int[][], Integer> findPeople() throws InterruptedException{
		int len = dataSize;
		V5FindPopulation[] partial = new V5FindPopulation[4];
		for (int i = 0; i < 4; i++){
      if(i != 3){
  			partial[i] = new V5FindPopulation(len / 4 * i, len / 4 * (i + 1), arr, dataSize, corners, cols, rows);
  			partial[i].start();
      } else {
        partial[i] = new V5FindPopulation(len / 4 * i, len, arr, dataSize, corners, cols, rows);
  			partial[i].start();
      }
		}

		for (int i = 0; i < 4; i++){
			partial[i].join();
		}
    //System.out.println(totalPopulation);
    /*for(int i = 0; i < lock.length; i++) {
      for(int j = 0; j < lock[0].length; j++) {
        System.out.println(i+ "," + j + ":" + matrixPeople[i][j]);
      }
    }*/
		return new Pair<int[][], Integer>(matrixPeople, totalPopulation);
	}
	@Override
	public void run(){

		for(int i = low; i < high; i++) {
      synchronized(lock){
        CensusGroup group = arr[i];
  			int[] index = findIndex(group.longitude, group.latitude, corners, cols, rows);
  			int rowNum = index[0];
  			int colNum = index[1];
  			if (rowNum == rows && colNum == cols){
  				matrixPeople[rowNum - 1][colNum - 1] += group.population;
  			}else if (rowNum == rows){
  				matrixPeople[rowNum - 1][colNum] += group.population;
  			}else if (colNum == cols){
  				matrixPeople[rowNum][colNum - 1] += group.population;
  			}else{
  				matrixPeople[rowNum][colNum] += group.population;
  			}
  			totalPopulation += group.population;
      }
    }
	}
	public int[] findIndex(float longitude, float latitude, float[] corners, int cols, int rows) {
		float top = corners[0];
		float bottom = corners[1];
		float left = corners[2];
		float right = corners[3];

		float height = (top - bottom)/rows;
		float width = (right - left)/cols;

		int[] index = new int[2];

		index[0] = (int)((latitude-bottom)/height);
		index[1] = (int)((longitude-left)/width);
		return index;
	}
}

/**
 * This class queries census data to find population densities in different
 * areas of the US.
 */
public class PopulationQuery{

	/**
	 * For parsing - the number of comma-separated fields on a given line.
	 */
	public static final int TOKENS_PER_LINE = 7;

	/**
	 * For parsing - zero-based index of the field containing the population of
	 * the current census group.
	 */
	public static final int POPULATION_INDEX = 4;

	/**
	 * For parsing - zero-based index of the field containing the latitude of
	 * the current census group.
	 */
	public static final int LATITUDE_INDEX = 5;

	/**
	 * For parsing - zero-based index of the field containing the longitude of
	 * the current census group.
	 */
	public static final int LONGITUDE_INDEX = 6;

	/**
	 * There should be only one fork/join pool per program, so this needs to be
	 * a static variable.
	 */
	private static ForkJoinPool fjPool = new ForkJoinPool();

	/**
	 * Array of census data parsed from the input file.
	 */
	private CensusData data;
	private int[][] matrixPeople;
	private int totalPopulation = 0;
	/**
	 * Initialize the query object by parsing the census data in the given file.
	 *
	 * @param filename
	 *            name of the census data file
	 */
	public PopulationQuery(String filename) {
		// Parse the data and store it in an array.
		this.data = parse(filename);
	}

	/**
	 * Parse the input file into a large array held in a CensusData object.
	 *
	 * @param filename
	 *            name of the file to be used as input.
	 * @return CensusData object containing the parsed data.
	 */
	private static CensusData parse(String filename) {
		CensusData result = new CensusData();

		try {
			BufferedReader fileIn = new BufferedReader(new FileReader(filename));

			/*
			 * Skip the first line of the file. After that, each line has 7
			 * comma-separated numbers (see constants above). We want to skip
			 * the first 4, the 5th is the population (an int) and the 6th and
			 * 7th are latitude and longitude (floats).
			 */

			try {
				/* Skip the first line. */
				String oneLine = fileIn.readLine();

				/*
				 * Read each subsequent line and add relevant data to a big
				 * array.
				 */
				while ((oneLine = fileIn.readLine()) != null) {
					String[] tokens = oneLine.split(",");
					if (tokens.length != TOKENS_PER_LINE)
						throw new NumberFormatException();
					int population = Integer.parseInt(tokens[POPULATION_INDEX]);
					result.add(population,
							Float.parseFloat(tokens[LATITUDE_INDEX]),
							Float.parseFloat(tokens[LONGITUDE_INDEX]));
				}
			} finally {
				fileIn.close();
			}
		} catch (IOException ioe) {
			System.err
					.println("Error opening/reading/writing input or output file.");
			System.exit(1);
		} catch (NumberFormatException nfe) {
			System.err.println(nfe.toString());
			System.err.println("Error in file format");
			System.exit(1);
		}
		return result;
	}

	/**
	 * Preprocess the census data for a run using the given parameters.
	 *
	 * @param cols
	 *            Number of columns in the map grid.
	 * @param rows
	 *            Number of rows in the map grid.
	 * @param versionNum
	 *            implementation to use
	 */
	public float[] findCorner1(){
		float top = Integer.MIN_VALUE;
		float bottom = Integer.MAX_VALUE;
		float left = Integer.MAX_VALUE;
		float right = Integer.MIN_VALUE;
		for (int i = 0; i < this.data.dataSize; i++){
			CensusGroup group = this.data.data[i];
			top = Math.max(top, group.latitude);
			bottom = Math.min(bottom, group.latitude);
			left = Math.min(left, group.longitude);
			right = Math.max(right, group.longitude);
		}

		return new float[]{top, bottom, left, right};
	}


	public void findTotalPeople1(int cols, int rows, float[] corners) {
		matrixPeople = new int[rows][cols];

		for(int i = 0; i < this.data.dataSize; i++) {
			CensusGroup group = this.data.data[i];
			int[] index = findIndex(group.longitude, group.latitude, corners, cols, rows);
			int rowNum = index[0];
			int colNum = index[1];
			if (rowNum == rows && colNum == cols){
				matrixPeople[rowNum - 1][colNum - 1] += group.population;
			}else if (rowNum == rows){
				matrixPeople[rowNum - 1][colNum] += group.population;
			}else if (colNum == cols){
				matrixPeople[rowNum][colNum - 1] += group.population;
			}else{
				matrixPeople[rowNum][colNum] += group.population;
			}
			totalPopulation += group.population;
		}

		//return matrixPeople;
	}
	public int[] findIndex(float longitude, float latitude, float[] corners, int cols, int rows) {
		float top = corners[0];
		float bottom = corners[1];
		float left = corners[2];
		float right = corners[3];

		float height = (top - bottom)/rows;
		float width = (right - left)/cols;

		int[] index = new int[2];
		index[0] = (int)((latitude-bottom)/height);
		index[1] = (int)((longitude-left)/width);
		return index;
	}

	public void preprocess(int cols, int rows, int versionNum) throws InterruptedException{
		// YOUR CODE GOES HERE
		float[] corners = new float[4];
		switch(versionNum){
			case 1:
				corners = findCorner1();
				findTotalPeople1(cols, rows, corners);
				break;
			case 2:
				//go here
				V2FindCorner v2FindCorner = new V2FindCorner(0, data.dataSize - 1, data.data, data.dataSize);
				corners = v2FindCorner.findCorner2();
				//(Arrays.toString(corners));
				V2FindPopulation v2FindPopulation = new V2FindPopulation(0, data.dataSize - 1, data.data, data.dataSize, corners, cols, rows);
				Pair<Integer, int[][]> pair = v2FindPopulation.findPopulation2();
				totalPopulation = pair.getElementA();
				matrixPeople = pair.getElementB();
				break;
			case 3:
				corners = findCorner1();
				findTotalPeople1(cols, rows, corners);
				v3preprocess();
			case 4:
				//for version 4, we don't need to write new codes, because we can reuse the
				//the code written before
				//version 4, parallelize the find corners method and the build grid method
				V2FindCorner v4FindCorner = new V2FindCorner(0, data.dataSize - 1, data.data, data.dataSize);
				corners = v4FindCorner.findCorner2();

				V2FindPopulation v4FindPopulation = new V2FindPopulation(0, data.dataSize - 1, data.data, data.dataSize, corners, cols, rows);
				Pair<Integer, int[][]> pair4 = v4FindPopulation.findPopulation2();
				totalPopulation = pair4.getElementA();
				//matrixpeople works as the distribution of the total population over grids
				matrixPeople = pair4.getElementB();

				//sequentially add up the population in step 2
				v3preprocess();
			case 5:
				V2FindCorner v5FindCorner = new V2FindCorner(0, data.dataSize, data.data, data.dataSize);
				corners = v5FindCorner.findCorner2();
        		V5FindPopulation.matrixPeople = new int[rows][cols];
				V5FindPopulation v5FindPopulation = new V5FindPopulation(0, data.dataSize, data.data, data.dataSize, corners, cols, rows);
				Pair<int[][], Integer> pair5 = v5FindPopulation.findPeople();
				totalPopulation = pair5.getElementB();
				//matrixpeople works as the distribution of the total population over grids
				matrixPeople = pair5.getElementA();
        
				v3preprocess();
			default:
				break;
		}
	}

	/**
	 * Query the population of a given rectangle.
	 *
	 * @param w
	 *            western edge of the rectangle
	 * @param s
	 *            southern edge of the rectangle
	 * @param e
	 *            eastern edge of the rectangle
	 * @param n
	 *            northern edge of the rectangle
	 * @return pair containing the population of the rectangle and the
	 *         population as a percentage of the total US population.
	 */
	public Pair<Integer, Float> singleInteraction(int w, int s, int e, int n) {
		// YOUR CODE GOES HERE
		int population = 0;
		float percentage = (float)0;
		for (int i = s; i <= n; i++){
			for (int j = w; j <= e; j++){
				population += matrixPeople[i - 1][j - 1];
			}
		}
		percentage = (float)population / (float)totalPopulation * 100;
		Pair<Integer, Float> pair = new Pair<Integer, Float>(population, percentage);
		return pair;

	}
	//VERSION 3: v3preprocess pre-compute the rectangle of population we need by applying the
	//dynamic programming algorithm shown on the site page
	public void v3preprocess(){
		//step1: init the first row and the first col
		int row = matrixPeople.length;
		int col = matrixPeople[0].length;

		for (int j = 1; j < col; j++){
			matrixPeople[row - 1][j] = matrixPeople[row - 1][j - 1] + matrixPeople[row - 1][j];
		}
		for (int i = row - 2; i >= 0; i--){
			matrixPeople[i][0] = matrixPeople[i + 1][0] + matrixPeople[i][0];
		}
		//step2: dp main algorithm
		for (int i = row - 2; i >= 0; i--){
			for (int j = 1; j < col; j++){
				matrixPeople[i][j] = matrixPeople[i + 1][j] + matrixPeople[i][j - 1] - matrixPeople[i + 1][j - 1] + matrixPeople[i][j];
			}
		}
	}
	//VERSION 3: we apply the O(1) query algorithm discussed on the site page
	public Pair<Integer, Float> singleInteractionV3(int w, int s, int e, int n){
		int row = matrixPeople.length;
		int col = matrixPeople[0].length;
		int bottomRight = matrixPeople[s-1][e-1];
		int topRight = (n == row) ? 0 : matrixPeople[n][e-1];
		int bottomLeft = (w == 1) ? 0 : matrixPeople[s-1][w-2];
		int topLeft = (n == row || w == 1) ? 0 : matrixPeople[n][w-2];
		int result = bottomRight - topRight - bottomLeft + topLeft;

		float percentage = (float)result / (float)totalPopulation * 100;
		return new Pair<Integer, Float>(result, percentage);
	}
	// argument 1: file name for input data: pass this to parse
	// argument 2: number of x-dimension buckets
	// argument 3: number of y-dimension buckets
	// argument 4: -v1, -v2, -v3, -v4, or -v5
	public static void main(String[] args) throws InterruptedException{
		// Parse the command-line arguments.
		String filename;
		int cols, rows, versionNum;
		try {
			filename = args[0];
			cols = Integer.parseInt(args[1]);
			rows = Integer.parseInt(args[2]);
			String versionStr = args[3];
			Pattern p = Pattern.compile("-v([12345])");
			Matcher m = p.matcher(versionStr);
			m.matches();
			versionNum = Integer.parseInt(m.group(1));
		} catch (Exception e) {
			System.out
					.println("Usage: java PopulationQuery <filename> <rows> <cols> -v<num>");
			System.exit(1);
			return;
		}

		// Parse the input data.
		PopulationQuery pq = new PopulationQuery(filename);

		// Preprocess the input data.
		pq.preprocess(cols, rows, versionNum);

		// Read queries from stdin.
		Scanner scanner = new Scanner(new BufferedInputStream(System.in));
		while (true) {
			int w, s, e, n;
			try {
				System.out.print("Query? (west south east north | quit) ");
				String west = scanner.next();
				if (west.equals("quit")) {
					break;
				}
				w = Integer.parseInt(west);
				s = scanner.nextInt();
				e = scanner.nextInt();
				n = scanner.nextInt();

				if (w < 1 || w > cols)
					throw new IllegalArgumentException();
				if (e < w || e > cols)
					throw new IllegalArgumentException();
				if (s < 1 || s > rows)
					throw new IllegalArgumentException();
				if (n < s || n > rows)
					throw new IllegalArgumentException();
			} catch (Exception ex) {
				System.out
						.println("Bad input. Please enter four integers separated by spaces.");
				System.out.println("1 <= west <= east <= " + cols);
				System.out.println("1 <= south <= north <= " + rows);
				continue;
			}
			Pair<Integer, Float> result;
			// Query the population for this rectangle.
			if(versionNum != 3 && versionNum != 4 && versionNum != 5){
				result = pq.singleInteraction(w, s, e, n);
			} else {
				//only call this part for version 3 and version 4
				//System.out.println("This is v");
				result = pq.singleInteractionV3(w, s, e, n);
			}
			System.out.printf("Query population: %10d\n", result.getElementA());
			System.out.printf("Percent of total: %10.2f%%\n",
					result.getElementB());
		}
		scanner.close();
	}
}
