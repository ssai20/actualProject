import java.util.function.BiFunction;
import java.util.function.Function;



class Derivative {
    int node;
    Function<Double, Double> function;
    Function<Double, Double> phi;
    Function<Double, Double> phiDerivative;
    Function<Double, Double> accuracyDerivative;
    BiFunction<Double, Double, Double> newMethodDerivative1;
    Function<Double, Double> newMethodDerivative2;
    public Derivative (int node, Function<Double, Double> function, Function<Double,Double> accuracyDerivative, Function<Double, Double> phi, Function<Double, Double> phiDerivative){
        this.node = node;
        this.function = function;
        this.accuracyDerivative = accuracyDerivative;
        this.phi = phi;
        this.phiDerivative = phiDerivative;
    }


    public double find (int N, double epsilon){
        int L = (node - 1) * (N - 1);
        double[] u = new double[L+1];
        double[] uAccuracyDerivative = new double[L+1];
        double[] Phi = new double[L+1];
        double[] PhiDerivative = new double[L+1];
        double[] uDerivativeNewMethod = new double[L+1];
        double h = 1./L;
        double[] x = new double[L+1];
        x[0] = 0.;
        for (int j = 1;j<=L; j++){
            x[j] = x[j-1] + h;
        }

        for (int i=0;i<=L;i++){
            u[i] = function.apply(x[i]);
        }
        for (int i=0;i<=L;i++) {
            uAccuracyDerivative[i] = accuracyDerivative.apply(x[i]);
        }
        for (int i=0;i<=L;i++) {
            Phi[i] = phi.apply(x[i]);
        }
        for (int i=0;i<=L;i++) {
            PhiDerivative[i] = phiDerivative.apply(x[i]);
        }
        for (int i=0;i<=L;i++) {
            uAccuracyDerivative[i] = accuracyDerivative.apply(x[i]);
        }

        for (int i = 0; i < L/(node-1); i++ ){
            for (int j = i * (node-1); j <= (node-1) * (i+1); j++){
                uDerivativeNewMethod[j] = (u[(2*i+1)*(node-1)/2]-u[i*(node-1)]) * PhiDerivative[j]/(Phi[(2*i+1)*(node-1)/2]-Phi[i*(node-1)]);
            }
        }

        double[] norm = new double[L+1];
        for (int i=0;i<=L;i++){
            norm[i] = epsilon*Math.abs(uDerivativeNewMethod[i] - uAccuracyDerivative[i]);
        }
        double maxNorm = 0.;
        for(int i=0;i<=L;i++){
            if (maxNorm<norm[i]) maxNorm = norm[i];
        }
        return maxNorm;
    }




}
public class ZX {

    public static double findUexp (int N, double epsilon) {
        int L = 4*(N-1);
        double[] U = new double[L+1];
        double[] proizvU = new double[L+1];
        double[] proizvPogranSloiU = new double[L+1];
        double hh = 1./L;
        double[] x = new double[L+1];

        x[0] = 0;
        for (int j=1;j<=L;j++){
            x[j] = x[j-1]+hh;
        }

        for (int i=0;i<=L;i++){
            U[i] = Math.cos(Math.PI*x[i]) + Math.exp(-x[i]/epsilon);
        }
        for (int i=1;i<=L;i++){
            proizvU[i] = (U[i] - U[i-1])/hh;
//            proizvPogranSloiU[i] = (-1./epsilon)*Math.exp(-x[i]/epsilon)*(U[i] - U[i-1])/
//                    (Math.exp(-x[i]/epsilon) - Math.exp(-x[i-1]/epsilon));
        }
        //вычисление первой производной на 5 точках
        for (int i=2;i<=L-2;i=i+4){
            proizvU[i-2] = (U[i] - U[i-1])/hh;
            proizvU[i-1] = (U[i] - U[i-1])/hh;
            proizvU[i] = (U[i] - U[i-1])/hh;
            proizvU[i+1] = (U[i] - U[i-1])/hh;
            proizvU[i+2] = (U[i] - U[i-1])/hh;
        }

        for (int i=2;i<=L-2;i=i+4){
            proizvPogranSloiU[i-2] = (-1./epsilon)*Math.exp(-x[i-2]/epsilon)*(U[i] - U[i-1])/
                    (Math.exp(-x[i]/epsilon) - Math.exp(-x[i-1]/epsilon));
            proizvPogranSloiU[i-1] = (-1./epsilon)*Math.exp(-x[i-1]/epsilon)*(U[i] - U[i-1])/
                    (Math.exp(-x[i]/epsilon) - Math.exp(-x[i-1]/epsilon));
            proizvPogranSloiU[i] = (-1./epsilon)*Math.exp(-x[i]/epsilon)*(U[i] - U[i-1])/
                    (Math.exp(-x[i]/epsilon) - Math.exp(-x[i-1]/epsilon));
            proizvPogranSloiU[i+1] = (-1./epsilon)*Math.exp(-x[i+1]/epsilon)*(U[i] - U[i-1])/
                    (Math.exp(-x[i]/epsilon) - Math.exp(-x[i-1]/epsilon));
            proizvPogranSloiU[i+2] = (-1./epsilon)*Math.exp(-x[i+2]/epsilon)*(U[i] - U[i-1])/
                    (Math.exp(-x[i]/epsilon) - Math.exp(-x[i-1]/epsilon));
        }

        double[] norm = new double[L+1];
/*
        for (int i=2;i<=L-2;i=i+4){
            norm[i-2] = Math.exp((x[i-1]-x[i-2])/epsilon) * (Math.cos(Math.PI*x[i]) - Math.cos(Math.PI*x[i-1]))/(1.-Math.exp(-hh/epsilon)) + epsilon*Math.PI*Math.sin(Math.PI*x[i-2]);
            norm[i-1] = Math.exp((x[i-1]-x[i-1])/epsilon) * (Math.cos(Math.PI*x[i]) - Math.cos(Math.PI*x[i-1]))/(1.-Math.exp(-hh/epsilon)) + epsilon*Math.PI*Math.sin(Math.PI*x[i-1]);
            norm[i] = Math.exp((x[i-1]-x[i])/epsilon) * (Math.cos(Math.PI*x[i]) - Math.cos(Math.PI*x[i-1]))/(1.-Math.exp(-hh/epsilon)) + epsilon*Math.PI*Math.sin(Math.PI*x[i]);
            norm[i+1] = Math.exp((x[i-1]-x[i+1])/epsilon) * (Math.cos(Math.PI*x[i]) - Math.cos(Math.PI*x[i-1]))/(1.-Math.exp(-hh/epsilon)) + epsilon*Math.PI*Math.sin(Math.PI*x[i+1]);
            norm[i+2] = Math.exp((x[i-1]-x[i+2])/epsilon) * (Math.cos(Math.PI*x[i]) - Math.cos(Math.PI*x[i-1]))/(1.-Math.exp(-hh/epsilon)) + epsilon*Math.PI*Math.sin(Math.PI*x[i+2]);
        }
*/

        double maxNorm = 0;

        for (int i=0;i<=L;i++){
//            norm[i] = epsilon*Math.abs(proizvU[i]-(-Math.PI*Math.sin(Math.PI*x[i]) - Math.exp(-x[i]/epsilon)/epsilon));

            norm[i] = epsilon*Math.abs(proizvPogranSloiU[i] -
                    (-Math.PI*Math.sin(Math.PI*x[i]) - Math.exp(-x[i]/epsilon)/epsilon));
        }

        for(int i=0;i<=L;i++){
            if (maxNorm<norm[i]) maxNorm = norm[i];
        }
        java.text.NumberFormat formatter = new java.text.DecimalFormat("0.##E0");
        String formattedDouble = formatter.format(maxNorm).replace(",", ".").replace("E", "e");
//            String formattedDouble = String.format("%6.2e", b).replace(",", ".");



        return maxNorm;
    }







    public static double findUsqrt (int N, double epsilon) {
        int L = 5*N;
        double[] U = new double[L+1];
        double[] proizvU = new double[L+1];
        double[] proizvPogranSloiU = new double[L+1];
        double hh = 1./L;
        double[] x = new double[L+1];

        System.out.println("first differentional x+epsilon:");
        x[0] = 0;
        for (int j=1;j<=L;j++){
            x[j] = x[j-1]+hh;
        }

        for (int i=0;i<=L;i++){
            U[i] = Math.cos(Math.PI*x[i]) + Math.exp(-x[i]/epsilon);
        }
        for (int i=1;i<=L;i++){
            proizvU[i] = (U[i] - U[i-1])/hh;
            proizvPogranSloiU[i] = (U[i] - U[i-1])/
                    ((Math.sqrt(x[i] + epsilon) - Math.sqrt(x[i-1] + epsilon))*(2.*Math.sqrt(x[i]+epsilon)));
        }
        double maxNorm = 0;
        double[] norm = new double[L+1];
        for (int i=1;i<=L;i++){
            norm[i] = Math.sqrt(epsilon)*Math.abs(proizvU[i]-(-Math.PI*Math.sin(Math.PI*x[i]) - Math.exp(-x[i]/epsilon)/epsilon));
//            norm[i] = epsilon*Math.abs(proizvPogranSloiU[i] -
//                    (-Math.PI*Math.sin(Math.PI*x[i]) - Math.exp(-x[i]/epsilon)/epsilon));
        }

        for(int i=1;i<=L;i++){
            if (maxNorm<norm[i]) maxNorm = norm[i];
        }
        return maxNorm;
    }




    public static double findSecondUexp (int N, double epsilon) {
        int L = 4*(N-1);
//        System.out.println(L);
        double[] U = new double[L+1];
        double[] Fi = new double[L+1];
        double[] secondProizvU = new double[L+1];
        double[] secondProizvPogranSloiU = new double[L+1];
//        double h = 1/N;
        double hh = 1./L;
        double[] uzelX = new double[N+1];
        double[] x = new double[L+1];

//        uzelX[0] = 0;
//        for (int i=1;i<=N;i++){
//            uzelX[i] = uzelX[i-1]+h;
//        }
        x[0] = 0;
        for (int j=1;j<=L;j++){
            x[j] = x[j-1]+hh;
//            System.out.println(x[j]);
        }

        for (int i=0;i<=L;i++){
            U[i] = Math.cos(Math.PI*x[i]) + Math.exp(-x[i]/epsilon);
            Fi[i] = Math.exp(-x[i]/epsilon);
//            System.out.println(U[i]);
        }
//        for (int i=1;i<=L-1;i=i+4){
//            secondProizvU[i-1] = (U[i+1] - U[i-1])/(hh*hh);
//            secondProizvU[i] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
//            secondProizvU[i+1] = (U[i+1] - U[i-1]) /(hh*hh);
//
////            System.out.println(proizvU[i]);
////            secondProizvPogranSloiU[i-1] = (1./(epsilon*epsilon))*Math.exp(-x[i]/epsilon)*(U[i+1] -2*U[i] + U[i-1])/(Math.exp(-x[i+1]/epsilon) -
////                    2*Math.exp(-x[i]/epsilon) + Math.exp(-x[i-1]/epsilon)); //false
//            secondProizvPogranSloiU[i] = (1./(epsilon*epsilon))*Math.exp(-x[i]/epsilon)*(U[i+1] - 2*U[i] + U[i-1])/(Math.exp(-x[i+1]/epsilon) -
//                    2*Math.exp(-x[i]/epsilon) + Math.exp(-x[i-1]/epsilon));
////            System.out.println(proizvPogranSloiU[i]);
//        }

/*

    //вычисление второй производной в 3 точках:
        for (int i=1;i<=L-1;i=i+2){
            secondProizvU[i-1] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i+1] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
        }
        for (int i=1;i<=L-1;i=i+2){
            secondProizvPogranSloiU[i-1] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i-1]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i+1] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i+1]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
        }
*/

        //вычисление второй производной на 5 точках
        for (int i=2;i<=L-2;i=i+4){
            secondProizvU[i-2] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i-1] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i+1] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i+2] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
        }
        for (int i=2;i<=L-2;i=i+4){
            secondProizvPogranSloiU[i-2] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i-2]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i-1] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i-1]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i+1] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i+1]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i+2] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i+2]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
        }

        double maxNorm = 0;
        double[] norm = new double[L+1];
        /*
        for (int i=0;i<=L;i++){
//            norm[i] = epsilon*epsilon*Math.abs(secondProizvPogranSloiU[i] - (-Math.PI*Math.PI*Math.cos(Math.PI*x[i]) +
//                    + Math.exp(-x[i]/epsilon)/(epsilon*epsilon)));
            norm[i] = epsilon*epsilon*Math.abs(secondProizvU[i] - (-Math.PI*Math.PI*Math.cos(Math.PI*x[i]) +
                    + Math.exp(-x[i]/epsilon)/(epsilon*epsilon)));

        }
        */

        for (int i=2;i<=L-2;i=i+4){
            norm[i-2] = Math.exp((x[i-1]-x[i-2])/epsilon) * (Math.cos(Math.PI*x[i-1])-2.*Math.cos(Math.PI*x[i])+Math.cos(Math.PI*x[i+1]))/(1.-2.*Math.exp(-hh/epsilon)+Math.exp(-2. * hh/epsilon)) +
                    epsilon*epsilon*Math.PI*Math.PI*Math.cos(Math.PI*x[i-2]);
            norm[i-1] = Math.exp((x[i-1]-x[i-1])/epsilon) * (Math.cos(Math.PI*x[i-1])-2.*Math.cos(Math.PI*x[i])+Math.cos(Math.PI*x[i+1]))/(1.-2.*Math.exp(-hh/epsilon)+Math.exp(-2. * hh/epsilon)) +
                    epsilon*epsilon*Math.PI*Math.PI*Math.cos(Math.PI*x[i-1]);
            norm[i] = Math.exp((x[i-1]-x[i])/epsilon) * (Math.cos(Math.PI*x[i-1])-2.*Math.cos(Math.PI*x[i])+Math.cos(Math.PI*x[i+1]))/(1.-2.*Math.exp(-hh/epsilon)+Math.exp(-2. * hh/epsilon)) +
                    epsilon*epsilon*Math.PI*Math.PI*Math.cos(Math.PI*x[i]);
            norm[i+1] = Math.exp((x[i-1]-x[i+1])/epsilon) * (Math.cos(Math.PI*x[i-1])-2.*Math.cos(Math.PI*x[i])+Math.cos(Math.PI*x[i+1]))/(1.-2.*Math.exp(-hh/epsilon)+Math.exp(-2. * hh/epsilon)) +
                    epsilon*epsilon*Math.PI*Math.PI*Math.cos(Math.PI*x[i+1]);
            norm[i+2] = Math.exp((x[i-1]-x[i+2])/epsilon) * (Math.cos(Math.PI*x[i-1])-2.*Math.cos(Math.PI*x[i])+Math.cos(Math.PI*x[i+1]))/(1.-2.*Math.exp(-hh/epsilon)+Math.exp(-2. * hh/epsilon)) +
                    epsilon*epsilon*Math.PI*Math.PI*Math.cos(Math.PI*x[i+2]);
        }

        for(int i=0;i<=L;i++){
            if (maxNorm<norm[i]) maxNorm = norm[i];
        }
        return maxNorm;
    }

    public static double findSecondUsqrt (int N, double epsilon) {
        int L = 3*N;
//        System.out.println(L);
        double[] U = new double[L+1];
        double[] Fi = new double[L+1];
        double[] secondProizvU = new double[L+1];
        double[] secondProizvPogranSloiU = new double[L+1];
//        double h = 1/N;
        double hh = 1./L;
        double[] uzelX = new double[N+1];
        double[] x = new double[L+1];

        x[0] = 0;
        for (int j=1;j<=L;j++){
            x[j] = x[j-1]+hh;
        }

        for (int i=0;i<=L;i++){
            U[i] = Math.cos(Math.PI*x[i]) + Math.exp(-x[i]/epsilon);
//            Fi[i] = Math.exp(-x[i]/epsilon);
            Fi[i] = Math.sqrt(x[i] + epsilon);
        }
        for (int i=1;i<=L-1;i=i+3){
            secondProizvU[i-1] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i+1] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
        }
        for (int i=1;i<=L-1;i=i+3){
//            secondProizvPogranSloiU[i] = (1./(epsilon*epsilon)) * Math.exp(-x[i]/epsilon)*(U[i+1] - 2*U[i] + U[i-1])/(Math.exp(-x[i+1]/epsilon) -
//                    2*Math.exp(-x[i]/epsilon) + Math.exp(-x[i-1]/epsilon));
            secondProizvPogranSloiU[i-1] = (U[i-1] - 2.*U[i] + U[i+1]) * (-1./(4.*(Math.sqrt((x[i-1]+epsilon)*(x[i-1]+epsilon)*(x[i-1]+epsilon)))))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i] = (U[i-1] - 2.*U[i] + U[i+1]) * (-1./(4.*(Math.sqrt((x[i]+epsilon)*(x[i]+epsilon)*(x[i]+epsilon)))))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i+1] = (U[i-1] - 2.*U[i] + U[i+1]) * (-1./(4.*(Math.sqrt((x[i+1]+epsilon)*(x[i+1]+epsilon)*(x[i+1]+epsilon)))))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
        }
        double maxNorm = 0;
        double[] norm = new double[L+1];
        for (int i=0;i<=L;i++){
            norm[i] = epsilon*epsilon*Math.abs(secondProizvPogranSloiU[i] - (-Math.PI*Math.PI*Math.cos(Math.PI*x[i]) +
                    + Math.exp(-x[i]/epsilon)/(epsilon*epsilon)));
//            norm[i] = epsilon*epsilon*Math.abs(secondProizvU[i] - (-Math.PI*Math.PI*Math.cos(Math.PI*x[i]) +
//                    + Math.exp(-x[i]/epsilon)/(epsilon*epsilon)));

        }

        for(int i=0;i<=L-1;i++){
            if (maxNorm<norm[i]) maxNorm = norm[i];
        }
        return maxNorm;
    }

    public static void latexInitial(){
        System.out.println("\\begin{table} [!htb]");
        System.out.println("\\caption { Погрешность классической формулы для вычисления второй производной в 5 точках}");
        System.out.println("\\begin{center}");
        System.out.println("\\begin{tabular}{c|c|c|c|c|c|c}");
        System.out.println("\\hline $\\varepsilon$ & \\multicolumn{6}{c}{$N$} \\\\");
        System.out.println("\\cline{2-7}& $32$&$64$& $128$&$256$&$512$&$1024$ \\\\");
    }
    public static void latexTable(String[][] residual, String[][] oa){
        System.out.println("\\hline $1$&$".concat(residual[0][0])+"$&$".concat(residual[0][1])+"$&$".concat(residual[0][2])+"$&$".concat(residual[0][3])+"$& $".concat(residual[0][4])+"$& $".concat(residual[0][5])+"$\\\\");
        System.out.println("$o.a.$ &".concat(oa[0][1])+"&".concat(oa[0][2])+"& ".concat(oa[0][3])+"&".concat(oa[0][4])+"&".concat(oa[0][5])+"&\\\\");
        System.out.println("$16^{-1}$&$".concat(residual[1][0])+"$&$".concat(residual[1][1])+"$&$".concat(residual[1][2])+"$&$".concat(residual[1][3])+"$&$".concat(residual[1][4])+"$& $".concat(residual[1][5])+"$\\\\");
        System.out.println("$o.a.$&".concat(oa[1][1])+"&".concat(oa[1][2])+"&".concat(oa[0][3])+"&".concat(oa[0][4])+"&".concat(oa[0][5])+"&\\\\");
        System.out.println("$32^{-1}$&$".concat(residual[2][0])+"$&$".concat(residual[2][1])+"$&$".concat(residual[2][2])+"$&$".concat(residual[2][3])+"$&$".concat(residual[2][4])+"$& $".concat(residual[2][5])+"$\\\\");
        System.out.println("$o.a.$&".concat(oa[2][1])+"&".concat(oa[2][2])+"&".concat(oa[2][3])+"&".concat(oa[2][4])+"&".concat(oa[2][5])+"&\\\\");
        System.out.println("$64^{-1}$&$".concat(residual[3][0])+"$&$".concat(residual[3][1])+"$&$".concat(residual[3][2])+"$&$".concat(residual[3][3])+"$&$".concat(residual[3][4])+"$& $".concat(residual[3][5])+"$\\\\");
        System.out.println("$o.a.$&".concat(oa[3][1])+"&".concat(oa[3][2])+"&".concat(oa[3][3])+"&".concat(oa[3][4])+"&".concat(oa[3][5])+"&\\\\");
        System.out.println("$128^{-1}$&$".concat(residual[4][0])+"$&$".concat(residual[4][1])+"$&$".concat(residual[4][2])+"$&$".concat(residual[4][3])+"$&$".concat(residual[4][4])+"$ & $".concat(residual[4][5])+"$\\\\");
        System.out.println("$o.a.$&".concat(oa[4][1])+"&".concat(oa[4][2])+"&".concat(oa[4][3])+"&".concat(oa[4][4])+"&".concat(oa[4][5])+"&  \\\\");
        System.out.println("$256^{-1}$&$".concat(residual[5][0])+"$&$".concat(residual[5][1])+"$&$".concat(residual[5][2])+"$&$".concat(residual[5][3])+"$&$".concat(residual[5][4])+"$& $".concat(residual[5][5])+"$\\\\");
        System.out.println("$o.a.$&".concat(oa[5][1])+"&".concat(oa[5][2])+"&".concat(oa[5][3])+"&".concat(oa[5][4])+"&".concat(oa[5][5])+"&\\\\");
        System.out.println("$512^{-1}$&$".concat(residual[6][0])+"$&$".concat(residual[6][1])+"$&$".concat(residual[6][2])+"$&$".concat(residual[6][3])+"$&$".concat(residual[6][4])+"$& $".concat(residual[6][5])+"$\\\\");
        System.out.println("$o.a.$&".concat(oa[6][1])+"&".concat(oa[6][2])+"&".concat(oa[6][3])+"&".concat(oa[6][4])+"&".concat(oa[6][5])+"&\\\\");


    }
    public static void latexEnd(){
        System.out.println("\\hline");
        System.out.println("        \\end{tabular}");
        System.out.println("    \\end{center}");
        System.out.println("\\end{table}");
    }

    public static void main(String[] args) {
        String residual[][] = new String[7][6];
        String oac[][] = new String[7][6];
//        double a = 0.;
        //double epsilon = 1./1.;
        int node = 5;
//        Derivative firstDerivative = new Derivative(node, x -> Math.cos(Math.PI * x) + Math.exp(-x/(epsilon)), x -> -Math.PI*Math.sin(Math.PI*x) - Math.exp(-x/epsilon)/epsilon, x -> Math.exp(-x/epsilon), x -> -Math.exp(-x/epsilon)/epsilon);
        latexInitial();
        int kRes=0;
        int kOa=0;
        for (int eps = 8;eps<=512;eps=eps*2) {
            if (eps==8.) eps = 1;
            double a = 0.;
            double epsilon = 1./eps;
//            Derivative firstDerivative = new Derivative(node, x -> Math.cos(Math.PI * x) + Math.exp(-x/(epsilon)), x -> -Math.PI*Math.sin(Math.PI*x) - Math.exp(-x/epsilon)/epsilon, x -> Math.exp(-x/epsilon), x -> -Math.exp(-x/epsilon)/epsilon);
            int lRes=0;
            int lOa=0;
            for (int i = 32; i <= 1024; i = 2 * i) {
                //            double b = findUexp(i, epsilon);
                Derivative firstDerivative = new Derivative(node, x -> Math.cos(Math.PI * x) + Math.exp(-x/(epsilon)), x -> -Math.PI*Math.sin(Math.PI*x) - Math.exp(-x/epsilon)/epsilon, x -> Math.exp(-x/epsilon), x -> -Math.exp(-x/epsilon)/epsilon);
                double b = firstDerivative.find(i, epsilon);
//                java.text.NumberFormat formatter = new java.text.DecimalFormat("0.##E0");
//                String formattedDouble = formatter.format(b).replace(",", ".").replace("E", "e");
                            String formattedDouble = String.format("%6.2e", b).replace(",", ".");
                double four = a / b;
                double oa = Math.log10(four) / Math.log10(2.);
                java.text.NumberFormat formatterOa = new java.text.DecimalFormat("0.#");
                String formattedDoubleOa = formatterOa.format(oa).replace(",", ".").replace("E", "e");
                a = b;
                residual[kRes][lRes] = formattedDouble;
                lRes++;
                oac[kOa][lOa] = formattedDoubleOa;
                lOa++;
//                System.out.println("epsilon = 1/" + 1. / epsilon);
//                System.out.println("Порядок точности log2 (||" + i / 2. + "||/||" + i + "||) = ".concat(formattedDoubleOa));
//                System.out.println("_______________________________________________________________");
//                System.out.println("||" + i + "|| = ".concat(formattedDouble));
            }
            kRes++;
            kOa++;
            if (eps==1.) eps = 8;
        }
        latexTable(residual, oac);
        latexEnd();
    }
}