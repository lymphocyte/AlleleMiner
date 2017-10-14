package main;

import java.io.IOException;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;

import algorithm.Compare_With_IMGT;
import tools.QualityController;

public class Filter {
    private static String vcf_mode = "download";
    private static String vcf_file_path = "";
    private static String gene_position_file_path = "";
    private static String prefix = "";
    private static String output_path = "";
    private static String database_path = "";
    private static boolean reversed;
    private static String genetype = "";
    private static String VJ_type = "";
    private static int threshold_num = 1;
    
    public static void main(String[] args) throws IOException, CompoundNotFoundException {
	// TODO Auto-generated method stub
	if (args.length != 0) {
	    ProcessArgs(args);
	    output_path = System.getProperty("user.dir") + "/" + prefix + "_alleles.txt";
	    database_path = System.getProperty("user.dir") + "/" + prefix + "_database.txt";
	} else {
	    ErrorMsg();
	    System.exit(0);
	}

//	 Processor processor = new Processor();
//	 System.out.println("path is " + output_path);
//	 processor.Process(vcf_file_path, gene_position_file_path,
//	 output_path, vcf_mode, threshold);
	/*
	 * Filter out the alleles and generate the new allele database
	 */
	QualityController Qc = new QualityController(output_path, database_path, false, threshold_num);

	Compare_With_IMGT cw = new Compare_With_IMGT(Qc.getAlleles(), genetype, VJ_type, 0);
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
		if (par.equals("-vcfmode")) {
		    vcf_mode = value;
		} else if (par.equals("-vcfpath")) {
		    vcf_file_path = value;

		} else if (par.equals("-geneposition")) {
		    gene_position_file_path = value;
		} else if (par.equals("-outputprefix")) {
		    prefix = value;
		} else if (par.equals("-threshold_num")) {
		    threshold_num = Integer.parseInt(value);
		} else if (par.equals("-VJ_type")) {
		    VJ_type = value;
		} else if (par.equals("-rev")) {
		    if (value.contains("yes") || value.contains("Yes")) {
			reversed = true;
		    } else {
			reversed = false;
		    }
		} else if (par.equals("-genetype")) {
		    genetype = value;
		} else {
		    ErrorMsg();
		    System.exit(0);
		}

	    }
	}

    }

    public static void ErrorMsg() {
	System.out.println("The use of AlleleMiner is: java -jar AlleleMiner.jar -parameter value");
	System.out.println(
		"Example use: java -jar AlleleMiner.jar -vcfmode local -vcfpath /Documents/file.vcf -geneposition /Documents/position.txt -outputprefix TRBV -threshold_num 1");
	System.out.println("-vcfmode. \"local\" or \"download\".");
	System.out.println("-vcfpath. The absolute path of the vcf file.");
	System.out.println("-geneposition. The position file of the genes.");
	System.out.println("-outputprefix. The name you want to give for your result file.");
	System.out.println(
		"-threshold_num. The minimum number of  the occurance of the haplotype to be regarded as an alelle, default is 1");

    }
    

}
