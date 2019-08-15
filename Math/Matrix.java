public class Matrix{
    
    private double[][] entries;

    //Dimensions of the matrix
    private int N, M;

    //Constructs a square matrix (identity matrix if 2nd arg is true, otherwise constructs a 0 matrix)
    public Matrix(int size, boolean identity){
        
        entries = new double[size][size];

        if(identity)
            for(int i = 0; i < size; i++)
                entries[i][i] = 1;
        
        N = entries.length;
        M = entries[0].length;
    }

    //Constructs a 0 matrix with specified dimensions
    public Matrix(int rows, int cols){
        entries = new double[rows][cols];
        N = entries.length;
        M = entries[0].length;
    }

    //Constructs a new matrix based on a preexisting 2D array (matrix is a deep copy)
    public Matrix(double[][] array){
        for(int i = 0; i < array.length; i++){
            if(array[i].length != array[0].length){
                System.out.println("MATRIX ERROR: Ensure 2D array used for matrix creation is rectangular");
                return;
            }
        }

        entries = new double[array.length][array[0].length];

        for(int i = 0; i < array.length; i++)
            System.arraycopy(array[i], 0, entries[i], 0, array[0].length);
            
        N = entries.length;
        M = entries[0].length;
    }

    //Constructs a column vector out of the given input
    public Matrix(double[] array){
        entries = new double[array.length][1];
        
        for(int i = 0; i < array.length; i++)
            entries[i][0] = array[i];

        N = entries.length;
        M = 1;
    }

    //Gets the number of rows this matrix has
    public int numRows(){
        return N;
    }

    //Gets the number of columns this matrix has
    public int numCols(){
        return M;
    }

    //Sets one element of the matrix to a value
    public void set(int row, int col, double value){
        entries[row][col] = value;
    }

    //Sets a row of elements to a given array
    public void setRow(int row, double[] array){
        if(array.length != M){
            System.out.println("MATRIX ERROR: Row replacement length mismatch");
            return;
        }

        System.arraycopy(array, 0, entries[row], 0, M);
    }

    //Sets a column of elements to a given array
    public void setCol(int col, double[] array){
        if(array.length != N){
            System.out.println("MATRIX ERROR: Column replacement length mismatch");
            return;
        }

        for(int i = 0; i < N; i++)
            entries[i][col] = array[i];

    }

    //Get the matrix element at specified row and column
    public double get(int row, int col){
        return entries[row][col];
    }

    //Get deep copy of a row
    public double[] getRow(int row){
        double[] r = new double[M];
        System.arraycopy(entries[row], 0, r, 0, M);
        return r;
    }

    //Get deep copy of a column
    public double[] getCol(int col){
        double[] c = new double[N];

        for(int i = 0; i < N; i++)
            c[i] = entries[i][col];
        
        return c;
    }

    //Returns deep copy of this array
    public Matrix copy(){
        return new Matrix(entries);
    }

    //Returns the determinant of this matrix
    public double det(){

        if(N != M){
            System.out.println("MATRIX ERROR: Cannot find determinant of Non-Square Matrix");
            return 0;
        }

        double[] row = new double[N];
        double total = 1;
        double determinant = 1;

        Matrix holder = this.copy();

        for(int i = 0; i < N; i++){
            int j = i;

            while(j < N && holder.get(j, i) == 0)
                j++;

            if(j == N)
                continue;

            if(j != i){

                for(int k = 0; k < N; k++)
                    holder.swap(j, k, i, k);
                
                byte sign = (byte) Math.pow(-1, j-i);
                determinant *= sign;
            }

            for(int k = 0; k < N; k++)
                row[k] = holder.get(i, k);

            for(int k = i + 1; k < N; k++){
                double diag = row[i];
                double nextRow = holder.get(k, i);

                for(int l = 0; l < N; l++)
                    holder.set(k, l, diag * holder.get(k, l) - nextRow * row[l]);

                total *= diag;
            }
        }

        for(int i = 0; i < N; i++)
            determinant *= holder.get(i, i);
            
        return determinant / total;
    }

    //Returns the transpose of this matrix
    public Matrix T(){
        double[][] transpose = new double[M][N];
        for(int i = 0; i < N; i++)
            for(int j = 0; j < M; j++)
                transpose[j][i] = entries[i][j];
        return new Matrix(transpose);
    }

    //Returns the sum of this matrix and another
    public Matrix add(Matrix other){
        
        if(N != other.numRows() || M != other.numCols()){
            System.out.println("MATRIX ERROR: Summand dimensions mismatched");
            return null;
        }

        double[][] sum = new double[N][M];
        for(int i = 0; i < sum.length; i++)
            for(int j = 0; j < sum[0].length; j++)
                sum[i][j] = entries[i][j] + other.get(i, j);

        return new Matrix(sum);
    }

    //Returns the difference of this matrix and another (other matrix is the negated one)
    public Matrix subtract(Matrix other){
        
        if(N != other.numRows() || M != other.numCols()){
            System.out.println("MATRIX ERROR: Subtractand dimensions mismatched");
            return null;
        }

        double[][] diff = new double[N][M];
        for(int i = 0; i < diff.length; i++)
            for(int j = 0; j < diff[0].length; j++)
                diff[i][j] = entries[i][j] - other.get(i, j);

        return new Matrix(diff);
    }

    //Returns another one of this matrix scaled by a factor
    public Matrix scale(double factor){
        double[][] scaled = new double[N][M];
        for(int i = 0; i < N; i++)
            for(int j = 0; j < M; j++)
                scaled[i][j] = entries[i][j] * factor;
        return new Matrix(scaled);
    }

    //Multiplies this matrix with another from the right; returns the product (standard algorithm)
    public Matrix mult(Matrix other){

        if(M != other.numRows()){
            System.out.println("MATRIX ERROR: Multiplication factors dimension mismatch");
            return null;
        }

        double[][] product = new double[N][other.numCols()];
        for(int i = 0; i < product.length; i++){
            for(int j = 0; j < product[0].length; j++){
                double dot = 0;
                for(int k = 0; k < other.numRows(); k++){
                    dot += entries[i][k] * other.get(k, j);
                }
                product[i][j] = dot;
            }
        }

        return new Matrix(product);
    }

    //Returns the inverse of this matrix (Gauss-Jordan Elimination)
    public Matrix invert(){

        if(N != M){
            System.out.println("MATRIX ERROR: Cannot invert a non-square matrix");
            return null;
        }

        if(det() == 0){
            System.out.println("MATRIX ERROR: Cannot invert a non-singular matrix");
            return null;
        }

        double temp = 0;
        Matrix inverse = copy();
        inverse.augmentIdentity();

        for(int i = inverse.numRows()-1; i > 0; i--)
            if(inverse.get(i-1, 0) < inverse.get(i, 0))
                inverse.swapRows(i, i-1);

        for(int i = 0; i < inverse.numRows(); i++){
            for(int j = 0; j < inverse.numRows(); j++){
                if(i != j){
                    temp = inverse.get(j, i) / inverse.get(i, i);

                    for(int k = 0; k < inverse.numCols(); k++)
                        inverse.set(j, k, inverse.get(j, k) - inverse.get(i, k) * temp);

                }
            }
        }

        for(int i = 0; i < inverse.numRows(); i++){
            temp = inverse.get(i, i);
            for(int j = 0; j < inverse.numCols(); j++){
                inverse.set(i, j, inverse.get(i, j) / temp);
            }
        }

        double[][] deAugmented = new double[N][N];
        for(int i = 0; i < N; i++)
            for(int j = M; j < inverse.numCols(); j++)
                deAugmented[i][j-M] = inverse.get(i, j);

        return new Matrix(deAugmented);
    }

    //Swaps placement of 2 elements
    public void swap(int row1, int col1, int row2, int col2){
        double temp = entries[row1][col1];
        entries[row1][col1] = entries[row2][col2];
        entries[row2][col2] = temp;
    }

    //Swaps placement of 2 specified rows
    public void swapRows(int row1, int row2){
        double[] temp = entries[row1];
        entries[row1] = entries[row2];
        entries[row2] = temp;
    }

    //Augments this matrix with the given matrix
    public void augment(Matrix m){
        if(N != m.numRows()){
            System.out.println("MATRIX ERROR: Augment matrix height mismatch");
            return;
        }

        for(int i = 0; i < N; i++){
            double[] newRow = new double[M + m.numCols()];
            System.arraycopy(entries[i], 0, newRow, 0, M);
            System.arraycopy(m.getRow(i), 0, newRow, M, m.numCols());
            entries[i] = newRow;
        }

        M += m.numCols();
    }

    //Augments this matrix with the identity
    public void augmentIdentity(){
        augment(new Matrix(N, true));
    }

    //Returns string representation of this matrix
    public String toString(){
        String matrixString = "{";
        for(int i = 0; i < N; i++){
            matrixString += "{";
            for(int j = 0; j < M; j++){
                matrixString += (entries[i][j] + (j == M - 1 ? "}" : ", "));
            }
            matrixString += i == N - 1 ? "}" : ",\n ";
        }
        return matrixString;
    }

}