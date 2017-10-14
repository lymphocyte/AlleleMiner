package main;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;

import algorithm.Compare_With_IMGT;
import tools.QualityController;

public class Run {
    private static String vcf_mode = "download";
    private static String vcf_file_path = "";
    private static String gene_position_file_path = "";
    private static String prefix = "";
    private static String output_path = "";
    private static String database_path = "";
    private static boolean reversed;
    private static String genetype = "";
    private static String VJ_type = "";
    private static String threshold = "0.000005";
    private static double threshold_num = 0.0003;

    public static void main(String[] args)
	    throws FileNotFoundException, IOException, SQLException, CompoundNotFoundException {
	// TODO Auto-generated method stub

	/*
	 * We generate 20 simulated database to test
	 */
	// for (int i = 1; i < 20; i++) {

	if (args.length != 0) {
	    ProcessArgs(args);
	    output_path = System.getProperty("user.dir") + "/" + prefix + "_alleles.txt";
	    database_path = System.getProperty("user.dir") + "/" + prefix + "_database.txt";
	} else {
	    ErrorMsg();
	    System.exit(0);
	}
//	
	Processor processor = new Processor();
//	System.out.println("path is " + output_path);
	processor.Process(vcf_file_path, gene_position_file_path, output_path, vcf_mode, threshold);
////	/*
////	 * Compare with IMGT test, Filter out the alleles and generate the new allele database.
////	 */
//	int threshold_numeric  = (int)(threshold_num * processor.GetTotalHaplotypes());
//	QualityController Qc = new QualityController(output_path, database_path, false, threshold_numeric);
//
//	Compare_With_IMGT cw = new Compare_With_IMGT(Qc.getAlleles(), genetype, VJ_type, 0);
	// }
    }

    public static void ProcessArgs(String[] arguments) {
	int length = arguments.length;
	// String parameters =
	// "-vcfmode-vcfpath-geneposition-outputprefix-threshold-rev-genetype-VJ_type";

	if (length % 2 != 0) {
	    ErrorMsg();
	    System.exit(0);
	} else {
	    for (int i = 0; i < length; i = i + 2) {
		String par = arguments[i];
		String value = arguments[i + 1];
		if (par.equals("-vcfmode") || par.equals("-m")) {
		    vcf_mode = value;
		} else if (par.equals("-vcfpath") || par.equals("-v")) {
		    vcf_file_path = value;

		} else if (par.equals("-geneposition") || par.equals("-p")) {
		    gene_position_file_path = value;
		} else if (par.equals("-outputprefix") || par.equals("-o")) {
		    prefix = value;
		} else if (par.equals("-threshold_num") || par.equals("-t")) {
		    threshold_num = Double.parseDouble(value);
		} else if (par.equals("-VJ_type") || par.equals("-e")) {
		    VJ_type = value;
		} else if (par.equals("-rev") || par.equals("-r")) {
		    if (value.contains("yes") || value.contains("Yes")) {
			reversed = true;
		    } else {
			reversed = false;
		    }
		} else if (par.equals("-genetype") || par.equals("-g")) {
		    genetype = value;
		} else {
		    ErrorMsg();
		    System.exit(0);
		}

	    }
	}

    }

    public static void ErrorMsg() {
	System.out.println("####################################################");
	System.out.println("####################################################");
	System.out.println("#####      	AlleleMiner         	     #######");
	System.out.println("#####      	Author: Yaxuan Yu      	     #######");
	System.out.println("#####      	Contact: yuyaxuan0@gmail.com #######");
	System.out.println("####################################################");
	System.out.println("####################################################");
	System.out.println(
		"--------------------------------------------------------------------------------------------------------------------------------------------------");
	System.out.println("The use of AlleleMiner is: java -jar AlleleMiner.jar -parameter value");
	System.out.println(
		"Example use: java -jar AlleleMiner.jar -vcfmode local -vcfpath /Documents/file.vcf -geneposition /Documents/position.txt -outputprefix TRBV -threshold_num 0.0003 -VJ_type V -genetype TRB");
	System.out.println(
		"--------------------------------------------------------------------------------------------------------------------------------------------------");

	System.out.println("-vcfmode(-m). \"local\" or \"download\".");
	System.out.println("-vcfpath(-v). The absolute path of the vcf file.");
	System.out.println("-geneposition(-p). The position file of the genes.");
	System.out.println("-outputprefix(-o). The name you want to give for your result file.");
	System.out.println("-VJ_type(-e). V or J genes? \"V\"or \"J\"");
	System.out.println("-rev(-r). if the gene position is reveresed? \"yes\" or \"no\"");
	System.out.println("-genetype(-g). The type of the gene, \"TRB\",\"IGH\",\"TRA\" and \"IGL\"");

	System.out.println(
		"-threshold_percentage(-t). The minimum proportion of  the occurance of the haplotype to be regarded as an alelle, default is 0.0003");
    }
}
