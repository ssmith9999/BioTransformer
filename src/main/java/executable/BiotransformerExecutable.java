package executable;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingArgumentException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.MissingOptionException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.btransformers.Cyp450BTransformer;
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.btransformers.EnvMicroBTransformer;
import biotransformer.btransformers.HGutBTransformer;
import biotransformer.btransformers.Phase2BTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.FileUtilities;
import biotransformer.utils.HumanSuperBioTransformer;
import biotransformer.utils.MetaboliteFinder;
import biotransformer.utils.MetaboliteFinder.FinderOption;
import biotransformer.utils.UniversalBioTransformer;
import biotransformer.utils.BiotransformerSequence;
import biotransformer.version.Version;

import org.apache.commons.io.FilenameUtils;


public class BiotransformerExecutable {

	
	
	public BiotransformerExecutable() {
		// TODO Auto-generated constructor stub
		
	}
	
	private static Options generateOptions(){
		
		final Option taskOption = Option.builder("k")
				.required(true)
				.hasArg(true)
				.argName("BioTransformer Task")
				.longOpt("task")
				.desc("The task to be permed: pred for prediction, or cid for compound identification ")
				.build();
		
		
		final Option biotransformerOption = Option.builder("b")
				.required(false)
				.hasArg(true)
				.argName("BioTransformer Option")
				.longOpt("btType")
				.desc("The type of description: Type of biotransformer - EC-based  (ecbased), CYP450 (cyp450), Phase II (phaseII), "
						+ "Human gut microbial (hgut), human super transformer* (superbio, or allHuman), Environmental microbial (envimicro)**.\n"
						+ "If option -m is enabled, the only valid biotransformer types are allHuman, superbio and env.")
				.build();

		final Option biotransformerSequenceOption = Option.builder("q")
				.required(false)
				.hasArg(true)
				.argName("BioTransformer Sequence Option")
				.longOpt("bsequence")
				.desc("Define an ordered sequence of biotransformer/nr_of_steps to apply. Choose only from the following "
						+ "BioTranformer Types: allHuman, cyp450, ecbased, env, hgut, and phaseII. For instance, the "
						+ "following string representation describes a sequence of 2 steps of CYP450 metabolism, followed by 1 "
						+ "step of Human Gut metabolism, 1 step of Phase II, and 1 step of Environmental Microbial Degradation:\n"
						+ "'cyp450:2; hgut:1; phaseII:1; env:1'")
				.build();

		
		final Option nrOfStepsOption = Option.builder("s")
				.required(false)
				.hasArg(true)
				.argName("Number of steps")
				.longOpt("nsteps")
				.desc("The number of steps for the prediction. This option can be set by the user for the EC-based, CYP450, Phase II, and Environmental microbial biotransformers. The default value is 1.")
				.build();
		
//		final Option inputOption = Option.builder("i")
//				.required(true)
//				.hasArg(true)
//				.argName("Input")
//				.longOpt("input")
//				.desc("The input, which can be a SMILES string, a Mol file, or SDF file.")
//				.build();
//		
//		final Option outputOption = Option.builder("o")
//				.required(false)
//				.hasArg(true)
//				.argName("Output")
//				.longOpt("ioutput")
//				.desc("The output file name (which must be a SDF file) for single input queries, or the output folder for multiple input queries.\n"
//						+ "When submitting a mutiple input query, an output folder must be provided.")
//				.build();


		final Option smiInputOption = Option.builder("ismi")
				.required(false)
				.hasArg(true)
				.argName("SMILES Input")
				.longOpt("ismiles")
				.desc("The input, which can be a SMILES string")
				.build();

		final Option molInputOption = Option.builder("imol")
				.required(false)
				.hasArg(true)
				.argName("MOL Input")
				.longOpt("molinput")
				.desc("The input, which can be a Mol file")
				.build();
		
		final Option sdfInputOption = Option.builder("isdf")
				.required(false)
				.hasArg(true)
				.argName("Sdf Input")
				.longOpt("sdfinput")
				.desc("The input, which can be an SDF file.")
				.build();
		
		final Option csvOutputOption = Option.builder("ocsv")
				.required(false)
				.hasArg(true)
				.argName("Csv Output")
				.longOpt("csvoutput")
				.desc("Select this option to return CSV output(s). You must enter an output filename")
				.build();
		
		
		final Option sdfOutputOption = Option.builder("osdf")
				.required(false)
				.hasArg(true)
				.argName("Sdf Output")
				.longOpt("sdfoutput")
				.desc("Select this option to return SDF output(s). You must enter an output filename")
				.build();
		
		

		final Option annotateOption = Option.builder("a")
				.required(false)
				.hasArg(false)
				.argName("Annotate")
				.longOpt("annotate")
				.desc("Search PuChem for each product, and store with CID and synonyms, when available.")
				.build();
		
		final Option indentificationMassOption = Option.builder("m")
				.required(false)
				.hasArg(true)
				.argName("Masses")
				.longOpt("masses")
				.desc("Semicolon-separated list of masses of compounds to identify")
				.build();
		
		final Option indentificationFormulaMetadataOption = Option.builder("f")
				.required(false)
				.hasArg(true)
				.argName("Formulas")
				.longOpt("formulas")
				.desc("Semicolon-separated list of formulas of compounds to identify")
				.build();
		
//		final Option indentificationMetadataOption = Option.builder("d")
//				.required(false)
//				.hasArg(true)
//				.argName("Metadata used for identification")
//				.longOpt("criteria")
//				.desc("The metadata that is used for identification. It must be wither 'mass', or 'formula', or 'combo' ")
//				.build();		
//		
//		final Option formulaOption = Option.builder("r")
//				.required(false)
//				.hasArg(true)
//				.argName("A semicolon-separated list of masses and/or formulas to identify")
//				.longOpt("dList")
//				.desc("Given the starting comound(s), Find all their metabolites with the specified masses, formulas, or combinations thereof, and show a metabolic pathway. "
//						+ "A semicolon-separated list of masses is expected. E.g. '308.08:C15H16O7;', 'C6H6O2;C15H16O7', '308.08:320.09'"
//						+ "If a combination of mass and formula is given, they shoulld be separated by a semicolon.")
//				.build();	
//		
//		final Option massformulaOption = Option.builder("l")
//				.required(false)
//				.hasArg(true)
//				.argName("List of masse/formulas combinations to identify")
//				.longOpt("mList")
//				.desc("Given the starting comound(s), Find all their metabolites with the specified masses and formulas, and show a metabolic pathway. A semicolon-separated list of masses is expected. E.g. '308.08:C15H16O7;C6H6O2;320.09'")
//				.build();			
		
		final Option massToleranceOption = Option.builder("t")
				.required(false)
				.hasArg(true)
				.argName("Mass Tolerance")
				.longOpt("mTolerance")
				.desc("Mass tolerance for metabolite identification (default is 0.01).")
				.build();

		final Option helpOption = Option.builder("h")
				.required(false)
				.hasArg(false)
				.argName("help")
				.longOpt("help")
				.desc("Prints the usage.")
				.build();
		
		final Options options = new Options();
		options.addOption(taskOption);
		options.addOption(biotransformerOption);
		options.addOption(biotransformerSequenceOption);
		options.addOption(nrOfStepsOption);
		options.addOption(smiInputOption);
		options.addOption(molInputOption);
		options.addOption(sdfInputOption);
		options.addOption(csvOutputOption);
		options.addOption(sdfOutputOption);
		options.addOption(indentificationMassOption);
		options.addOption(indentificationFormulaMetadataOption);	
		options.addOption(massToleranceOption);
		options.addOption(annotateOption);
		options.addOption(helpOption);

		return options;
	}
	
	public static CommandLine generateCommandLine(
			final Options options, final String[] commandLineArguments) throws ParseException{
		final CommandLineParser cmdLineParser = new DefaultParser();
		CommandLine commandLine = null;
		
//		System.out.println(Version.current);
		String header = "\nThis is the version " + Version.current + " of BioTransformer. BioTransformer is a software tool that predicts small molecule metabolism in mammals, their gut microbiota,"
				+ " as well as the soil/aquatic microbiota. BioTransformer also assists scientists in metabolite identification, based on the metabolism prediction. \n\n";
		
		String footer = "\n(* ) While the 'superbio' option runs a set number of transformation steps in a pre-defined order (e.g. deconjugation first, then Oxidation/reduction, etc.),"
				+ " the 'allHuman' option predicts all possible metabolites from any applicable reaction(Oxidation, reduction, (de-)conjugation) at each step."
				+ "(** ) For the environmental microbial biodegradation, all reactions (aerobic and anaerobic) are reported, and not only the aerobic biotransformations (as per default in the EAWAG BBD/PPS system)."
				+ "\n\n*********\n"
				+ "Examples:\n"
				+ "*********\n\n"
				+"1) To predict the biotransformation of a molecule from an SDF input using the human super transformer (option superbio) and annotate the metabolites with names and database IDs (from PubChem), run\n"
				+ "\n	java -jar biotransformer-" + Version.current +".jar -k pred -b superbio -isdf #{input file name} -osdf #{output file} -a."
				+ "\n\n2) To predict the 2-step biotransformation of Thymol (a monoterpene) using the human super transformer (option allHuman) using the SMILES input and saving to a CSV file, run"
				+ "\n	java -jar biotransformer-" + Version.current +".jar  -k pred -b allHuman -ismi \"CC(C)C1=CC=C(C)C=C1O\" -ocsv #{replace with output file name} -s 2"
				+ "\n\n3) Identify all human metabolites (max depth = 2) of Epicatechin (\"O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1\") with masses 292.0946 Da and 304.0946 Da, with a mass tolerance of 0.01 Da."
				+ " Provide an annotation (Common name, synonyms, and PubChem CID), when available."
				+ "\n	java -jar biotransformer-" + Version.current +".jar  -k cid -b allHuman -ismi \"O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1\" -osdf #{replace with output file name} -s 2 -m \"292.0946;304.0946\" -t 0.01 -a"
				+ "\n	- DO NOT forget the quotes around the SMILES string or the list of masses"
				+ "\n"
				+ "4) Simulate an order sequence of metabolism of Atrazine (\"CCNC1=NC(=NC(=N1)Cl)NC(C)C\"), starting with two steps of Cyp450 oxidation, followed by one step of conjugation."
				+ "\n java -jar biotransformer-1.1.4.jar -ismi \"CCNC1=NC(=NC(=N1)Cl)NC(C)C\" -osdf ~/atrazine-sequence.sdf -k pred -q \"cyp450:2; phaseII:1\"\n"
				+ "\nTo report issues, provide feedback, or ask questions, please send an e-mail the following address: djoumbou@ualberta.ca\n\n"
				+ "BioTransformer is offered to the public as a freely acessible software package under the GNU License GPL v2.1.Users are free"
				+ " to copy and redistribute the material in any medium or format. Moreover, they could modify, and build upon the material unfer "
				+ "the condition that they must give appropriate credit, provide links to the license, and indicate if changes were made. Furthermore, "
				+ "the above copyright notice and this permission notice must be included. Use and re-distribution of the these resources, in whole or in part, "
				+ "for commercial purposes requires explicit permission of the authors. We ask that all users of the BioTransformer software tool, the BioTransformer web server, "
				+ "or BioTransformerDB to cite the BioTransformer reference in any resulting publications, and to acknowledge the authors."
				+ "\n\n"
				+ "Important Notice:\n"
				+ "-----------------\n"
				+ "\t\tBioTransformer's environmental microbial degradation module uses data from the EAWAG's Biodegradation and Biocatalysis Database, which is licensed by EnviPath (https://envipath.com/license/)"
				+ "under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0). Therefore users of the environmental microbial degradation "
				+ "module must cite the following paper:\n"
				+ "\tenviPath–The environmental contaminant biotransformation pathway resource. J Wicker, T Lorsbach, M Gütlein, E Schmid, D Latino, S Kramer, K Fenner. Nucleic Acids Research, gkv1229 (full text)."
				+ "\n\nTo use the environmental microbial module for commercial purposes, users must request an apropriate license from EnviPath.";
		
//		String footer = "* While the superbio option runs a set number of transformation steps in a pre-defined order (e.g. deconjugation first, then Oxidation/reduction, etc.), "
//				+ "the allHuman option predicts all possible metabolites from any applicable reaction(Oxidation, reduction, (de-)conjudation) at each step."
//				+ "\n\n** For the environmental microbial biodegradation, all reactions (aerobic and anaerobic) are reported, and not only the aerobic biotransformations "
//				+ "(as per default in the EAWAG BBD/PPS system)."
//				+ "\n\n*********\nExamples:\n*********\n\n1) To predict the biotransformation of a molecule from an SDF input using the human super transformer, use java -jar biotransformer-" + Version.current +".jar -k pred -b superbio "
//				+ "-f sdf -i #{input file name} -o #{output folder}.\n\n"
//				+ "\t\t- For each of the query molecule in the input file, an outputfile will be created with the list of corresponding metabolites.\n\n"
//				+ "2) To predict the 2-step biotransformation of Thymol (a monoterpene) using the human super transformer (option allHuman) using the SMILES input, run\n\n"
//				+ "java -jar biotransformer-" + Version.current +".jar -k pred -b allHuman -f smiles -i \"CC(C)C1=CC=C(C)C=C1O\" -o #{replace with output file name} -s 2\n\n"
//				+ "Currently, the outputfile is SDF per default.\n\n"
//				+ "3) Identify all human metabolites (max depth = 2) of Epicatechin (\"O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1\") with masses 292.0946 Da and 304.0946 Da, with a mass tolerance "
//				+ "of 0.01 Da. Provide an annotation (Common name).\n\n"
//				+ "java -jar biotransformer-" + Version.current +".jar -k cid -b allHuman -f smiles -i \"O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1\" -o #{replace with output file name} "
//				+ "-s 2 -m \"292.0946; 304.0946\" -t 0.01 -a\n\t\t- DO NOT forget the quotes around the SMILES string or the list of masses.\n\n"
//				+ "To report issues, provide feedback, or ask questions, please send an e-mail the following address: djoumbou@ualberta.ca\n\n"
//				+ "BioTransformer is offered to the public as a freely acessible software package. Beside the prediction software, a manually curated database called BioTransformerDB is also available. The package is available under the GNU license GPL v2.1\n\n"
//				+ "Users are free to copy and redistribute the material in any medium or format. Moreover, they could modify, and build upon the material under the condition that they must give appropriate "
//				+ "credit, provide links to the license, and indicate if changes were made. Furthermore, use and re-distribution of the these resources, in whole or in part, for commercial purposes requires "
//				+ "explicit permission of the authors. We ask that all users of the BioTransformer software tool or BioTransformerDB to cite the BioTransformer reference in any resulting publications, and to acknowledge the authors.\n\n"
//				+ "Cite: Djoumbou-Feunang Y, Fiamoncini J, de la Fuente AG, Manach C, Greiner R, and Wishart DS; BioTransformer: A Comprehensive Computational Tool for Small Molecule Metabolism Prediction and Metabolite Identification; Journal of Cheminformatics 201911:2; DOI: 10.1186/s13321-018-0324-5.";
		
		HelpFormatter formatter = new HelpFormatter();

		
		try{
			commandLine = cmdLineParser.parse(options, commandLineArguments);
		}
		catch (MissingOptionException missingOptionException){
			
			if( Arrays.asList(commandLineArguments).contains("-h") || Arrays.asList(commandLineArguments).contains("--help")){
				formatter.printHelp("\njava -jar biotransformer-" + Version.current +".jar", header, options, footer, true);
			}
			else {
				System.out.println(missingOptionException.getLocalizedMessage());
			}			
		}
		catch (ParseException parseException){
			System.out.println("Could not parse the command line arguments "
					+ Arrays.toString(commandLineArguments) + "\nfor the following reaons:" 
					+ parseException);		
//			throw (parseException);
		}
		return commandLine;
	}
	
	
	
	public static void main(String[] args) throws Exception{
		IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
		SmilesParser	smiParser		= new SmilesParser(builder);

		Options options = generateOptions();
		CommandLine commandLine = generateCommandLine(options, args);

//		System.out.println(commandLine.getOptionValue("b"));
//		System.out.println(commandLine.getOptionValue("f").length());
//		System.out.println(commandLine.getOptionValue("m"));
		
		IAtomContainer singleInput = null;
		String inputFileName = null;
//		IAtomContainerSet containers = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		int nrOfSteps = 1;
		
		boolean annotate = false;
		String masses = null;
		String formulas = null;
		String metadata_input = null;
		Double defaultMassToleranceThreshold = 0.01;
		Double massToleranceThreshold = null;
		String task = null;
		FinderOption opt = null;
		String iFormat = null;
		String oFormat = null;
		String outputF = null;
		BiotransformerSequence biotransformerSeqeuence = null;
		double scoreThreshold = 0.5;
		
		int number_of_molecules = 0;
		int successful_predictions = 0;


		
		if(Arrays.asList(args).contains("-a") || Arrays.asList(args)
				.contains("--annotate")){		
			annotate = true;
			
			System.out.println("\n\n=============================================================");
			System.out.println("Compounds will be annotated using the PubChem API. Make sure");
			System.out.println("to have a secure internet connection.");
			System.out.println("=============================================================\n");
		}
		
		if(Arrays.asList(args).contains("-ismi") || Arrays.asList(args)
				.contains("--ismiles")){
			iFormat = "smi";
		}
		else if(Arrays.asList(args).contains("-imol") || Arrays.asList(args)
				.contains("--sdfinput")){
			iFormat = "mol";
		}
		else if(Arrays.asList(args).contains("-isdf") || Arrays.asList(args)
				.contains("--sdfinput")){
			iFormat = "sdf";
		}
				
		if(Arrays.asList(args).contains("-ocsv") || Arrays.asList(args)
				.contains("--csvoutput")){
			oFormat = "csv";
			outputF = commandLine.getOptionValue("ocsv").trim();
		}
		else if(Arrays.asList(args).contains("-osdf") || Arrays.asList(args)
				.contains("sdfoutput")){
			oFormat = "sdf";
			outputF = commandLine.getOptionValue("osdf").trim();
		}		
		
		if(commandLine !=null){
						
			if(commandLine.getOptionValue("k") != null){
				task = commandLine.getOptionValue("k").trim();
				if( !(task.contentEquals("pred") || task.contentEquals("cid")) ){
					throw new IllegalArgumentException("Invalid task(\"" +  task + "\") entered. Enter either 'pred' (prediction) or 'cid' (compound identification)");
				}				
			}
			else {
				throw new IllegalArgumentException("\n\tThe task type is missing. You must select either 'pred' (prediction) or 'cid' (compound identification)");
			}

			if(commandLine.getOptionValue("m") != null){
				masses = commandLine.getOptionValue("m").trim();
				if(masses.length() == 0){
					throw new MissingArgumentException("\n\tPlease enter a list of monoisotopic masses.");
				}			
//				System.out.println("MASS: " + masses);
			}
			
			if(commandLine.getOptionValue("f") != null){
				formulas = commandLine.getOptionValue("f").trim();
				if(formulas.length() == 0){
					throw new MissingArgumentException("\n\tPlease enter a list of chemical formulas.");
				}			
//				System.out.println("MASS: " + masses);
			}
			
//			System.out.println("ANNOTATE: " + annotate);
			if(task.contentEquals("cid")){
				if(commandLine.getOptionValue("t") != null){
					if(commandLine.getOptionValue("t").trim().length() == 0){
						throw new MissingArgumentException("\n\tThe option '-t' was used but the mass tolerance threshold is missing. Add a value or omit '-t' to use the default value (" + defaultMassToleranceThreshold +")");
					}
					else{
						massToleranceThreshold = Double.valueOf(commandLine.getOptionValue("t").trim());	
//						System.out.println("massTolerance: " + massTolerance);
					}
				}
					
				if(masses == null && formulas == null){
					throw new IllegalArgumentException("\n\tIdentification metadata are missing. Please add a list of masses (-m) or a list of formulas (-r)");
				}
				else if(formulas != null) {
					if(masses != null) {
						throw new IllegalArgumentException("\tList of masses and formulas provided. Please exclusively provide either a list of masses (-m) or a list of formulas (-r)");											
					}
					else if(massToleranceThreshold != null) {
						throw new IllegalArgumentException("\n\tA mass tolerance threshold is accepted only for mass-based, and not formula-based identification tasks. Please remove this argument.");											
						
					}
					opt = FinderOption.FORMULA;
					metadata_input = formulas;
//					System.out.println(opt + "\t" + metadata_input);
					/**
					 * Here, there user has not provided a mass threshold tolerance, but we set it to the default, just to avoid the NullPointerException
					 * The massToleranceThreshold argument will not be used for the FORMULA option.
					 */
					massToleranceThreshold = defaultMassToleranceThreshold;
				}
				else if(masses != null){
					opt = FinderOption.MASS;
					metadata_input = masses;
					if(massToleranceThreshold == null) {
						massToleranceThreshold = defaultMassToleranceThreshold;
					}
				}
			}

			final String biotransformerType = commandLine.getOptionValue("b");
			
			if(commandLine.getOptionValue("s") != null){
				nrOfSteps = Integer.valueOf(commandLine.getOptionValue("s"));
//				System.out.println("nrOfSteps: " + nrOfSteps);
			}	
			
			String bseq = commandLine.getOptionValue("q");
			if(bseq != null) {
				if(bseq.contentEquals("")) {
					throw new IllegalArgumentException("\n\tIllegalArgumentException. Please provide a valid biotranformer sequence.");				
				}
				else if(bseq.length()>1) {
					biotransformerSeqeuence = new BiotransformerSequence(bseq);
				}		
			}
			
			if(biotransformerType != null && commandLine.getOptionValue("q") != null) {
				throw new IllegalArgumentException("IllegalArgumentException: The parameters '-b' and '-seq' are mutually "
						+ "exclusive. While '-b' describes a specific biotransformer type, '-q' describes an order sequence of biotransformers to be applied.");
			}
			else if(biotransformerType == null && commandLine.getOptionValue("q") == null){
				throw new MissingArgumentException("MissingArgumentException: Please select either '-b' for a bitransformer type, or '-seq' "
						+ "for a biotransformer sequence.");
			}
			if(commandLine.getOptionValue("s") != null && commandLine.getOptionValue("q") != null) {
					throw new IllegalArgumentException("IllegalArgumentException: The parameters '-s' and '-seq' are mutually "
							+ "exclusive. '-q' describes an order sequence of biotransformers to be applied, each for a specific number of steps");
			}					

			 
//			final String format = commandLine.getOptionValue("f").toString().trim();
//			final String outputF = commandLine.getOptionValue("o");
			
			if(iFormat.contentEquals("smi")){
				String smi = commandLine.getOptionValue("ismi");
				singleInput = smiParser.parseSmiles(smi);			
			}
			else if(iFormat.contentEquals("mol")){
				inputFileName = commandLine.getOptionValue("imol");
//				containers = FileUtils.parseSdf(inputFileName);
				if(inputFileName == null){
					throw new MissingOptionException("\n\tPlease specify an input file name (Molfile or SDF). For more information, type java -jar biotransformer-" + Version.current +".jar --help.");
				}
//				if(outputF == null){
//					throw new MissingOptionException("A destination folder must be specified when your query molecules are provided in a file (Molfile or SDF). For more information, type java -jar biotransformer-" + Version.current +".jar --help.");
//				}
			}
			
			else if(iFormat.contentEquals("sdf")){
				inputFileName = commandLine.getOptionValue("isdf");
//				containers = FileUtils.parseSdf(inputFileName);
				if(inputFileName == null){
					throw new MissingOptionException("\n\tPlease specify an input file name (Molfile or SDF). For more information, type java -jar biotransformer-" + Version.current +".jar --help.");
				}
			}
			else {
				throw new IllegalArgumentException("\n\tInvalid input format option(" + iFormat + ") entered. It must be one of 'ismi','imol', or 'isdf'. Type java -jar biotransformer-" + Version.current +".jar --help for help.");
			}

			if(oFormat.contentEquals("csv")){
				outputF = commandLine.getOptionValue("ocsv");
				if(outputF == null){
					throw new MissingOptionException("\n\tPlease specify an output file name (CSV). For more information, type java -jar biotransformer-" + Version.current +".jar --help.");
				}
			}			
			else if(oFormat.contentEquals("sdf")){
				outputF = commandLine.getOptionValue("osdf");
				if(outputF == null){
					throw new MissingOptionException("\n\tPlease specify an output file name (SDF). For more information, type java -jar biotransformer-" + Version.current +".jar --help.");
				}
			}
			else {
				throw new IllegalArgumentException("\n\tInvalid output format option(" + oFormat + ") entered. It must be one of 'ocsv' or 'osdf'. Type java -jar biotransformer-" + Version.current +".jar --help for help.");
			}
			
//			System.out.println("TASK: "+ task.contentEquals("cid"));
//			System.out.println("TASK: "+ task.trim().contentEquals("cid"));
			
			if(task.contentEquals("cid")){
				
				if(metadata_input != null){
//					System.out.println("IDENTIFICATION TASK");
					
					if(biotransformerType != null) {
						if(biotransformerType.contentEquals("allHuman")){
							
							String[] mArr = metadata_input.trim().split(";");
							ArrayList<String> dmassesOrFormulas = new ArrayList<String>();
							
							for(int k = 0; k < mArr.length; k++){
								try{
									dmassesOrFormulas.add(mArr[k].trim());
								}
								catch(Exception e){
									System.err.println(e.getMessage());
								}
							}
							
							
							if (singleInput !=null){
								MetaboliteFinder mtf = new MetaboliteFinder();
								number_of_molecules++;
								System.out.println("\n\nMolecule no. " + number_of_molecules);
								if(oFormat.contentEquals("csv")){
									mtf.findAllHumanMetabolitesToCSV(singleInput, dmassesOrFormulas, massToleranceThreshold, nrOfSteps, annotate, outputF, opt);
								}
								else if(oFormat.contentEquals("sdf")){
									mtf.findAllHumanMetabolites(singleInput, dmassesOrFormulas, massToleranceThreshold, nrOfSteps, annotate, outputF, opt);
								}
							}
							else {
								MetaboliteFinder mtf = new MetaboliteFinder();
								IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
								IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
	
								for(IAtomContainer atc : containers.atomContainers()){
									number_of_molecules++;
									System.out.println("\n\nMolecule no. " + number_of_molecules);
									try {
										metabolites.add(mtf.findAllHumanMetabolites(atc, dmassesOrFormulas, massToleranceThreshold, nrOfSteps, annotate, opt));
									}
									catch(Exception e) {
										System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getLocalizedMessage());
									}
									
								}
								
								if(oFormat.contentEquals("csv")){
									FileUtilities.saveAtomContainerSetToCSV(metabolites, outputF);
								}
								else if(oFormat.contentEquals("sdf")){							
									SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputF));		
									sdfWriter.write(metabolites);
									sdfWriter.close();
								}
							}	
							
						}
						else if(biotransformerType.contentEquals("env")){
							String[] mArr = metadata_input.trim().split(";");
							ArrayList<String> dmassesOrFormulas = new ArrayList<String>();
							
							for(int k = 0; k < mArr.length; k++){
								try{
									dmassesOrFormulas.add(mArr[k].trim());
								}
								catch(Exception e){
									System.err.println(e.getMessage());
								}
							}
							
							
							if (singleInput !=null){
								number_of_molecules++;
								System.out.println("\n\nMolecule no. " + number_of_molecules);
								MetaboliteFinder mtf = new MetaboliteFinder();
								
								if(oFormat.contentEquals("csv")){
									mtf.findAllEnvMicroMetabolitesToCSV(singleInput, dmassesOrFormulas, massToleranceThreshold, nrOfSteps, annotate, outputF, opt);
								}
								else if(oFormat.contentEquals("sdf")){
									mtf.findAllEnvMicroMetabolites(singleInput, dmassesOrFormulas, massToleranceThreshold, nrOfSteps, annotate, outputF, opt);
								}
							}
							else {
								MetaboliteFinder mtf = new MetaboliteFinder();
								IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
								IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
	
								for(IAtomContainer atc : containers.atomContainers()){
									number_of_molecules++;
									System.out.println("\n\nMolecule no. " + number_of_molecules);
									try {
										metabolites.add(mtf.findAllEnvMicroMetabolites(atc, dmassesOrFormulas, massToleranceThreshold, nrOfSteps, annotate, opt));
									}
									catch(Exception e) {
										System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getLocalizedMessage());
									}
									
								}							
								if(oFormat.contentEquals("csv")){
									FileUtilities.saveAtomContainerSetToCSV(metabolites, outputF);
								}
								else if(oFormat.contentEquals("sdf")){							
									SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputF));		
									sdfWriter.write(metabolites);
									sdfWriter.close();
								}
							}					
						}
	
						else if(biotransformerType.contentEquals("superbio")){
	//						System.out.println(opt + "\t\t\t" + metadata_input);
							String[] mArr = metadata_input.trim().split(";");
							ArrayList<String> dmassesOrFormulas = new ArrayList<String>();
							
							for(int k = 0; k < mArr.length; k++){
								try{
									dmassesOrFormulas.add(mArr[k].trim());
								}
								catch(Exception e){
									System.err.println(e.getMessage());
								}
							}
								
							if (singleInput !=null){
								number_of_molecules++;
								System.out.println("\n\nMolecule no. " + number_of_molecules);
								MetaboliteFinder mtf = new MetaboliteFinder();
								
								if(oFormat.contentEquals("csv")){
									mtf.findSuperbioMetabolitesToCSV(singleInput, dmassesOrFormulas, massToleranceThreshold, annotate, outputF, opt);
								}
								else if(oFormat.contentEquals("sdf")){
									mtf.findSuperbioMetabolites(singleInput, dmassesOrFormulas, massToleranceThreshold, annotate, outputF, opt);
								}
								
								
							}
							else {
								MetaboliteFinder mtf = new MetaboliteFinder();
								IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
								IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
	
								for(IAtomContainer atc : containers.atomContainers()){
									number_of_molecules++;
									System.out.println("\n\nMolecule no. " + number_of_molecules);
									try {
										metabolites.add(mtf.findSuperbioMetabolites(atc, dmassesOrFormulas, massToleranceThreshold, annotate, opt));
									}
									catch(Exception e) {
										System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getLocalizedMessage());
									}
									
								}							
								if(oFormat.contentEquals("csv")){
									FileUtilities.saveAtomContainerSetToCSV(metabolites, outputF);
								}
								else if(oFormat.contentEquals("sdf")){							
									SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputF));		
									sdfWriter.write(metabolites);
									sdfWriter.close();
								}
	
							}					
						}
						else{
							throw new IllegalArgumentException("\n\tFor metabolite identification, the biotransformer type must be either allHuman, superbio, or env.");
						}
					}
					else if(biotransformerSeqeuence != null) {
						String[] mArr = metadata_input.trim().split(";");
						ArrayList<String> dmassesOrFormulas = new ArrayList<String>();
						
						for(int k = 0; k < mArr.length; k++){
							try{
								dmassesOrFormulas.add(mArr[k].trim());
							}
							catch(Exception e){
								System.err.println(e.getMessage());
							}
						}		
						IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
						MetaboliteFinder mtf = new MetaboliteFinder();
						if (singleInput !=null){
							number_of_molecules++;
							
							metabolites.add(mtf.findMetabolitesFromSequence(singleInput, 
									biotransformerSeqeuence, dmassesOrFormulas, massToleranceThreshold, annotate, opt));

							successful_predictions++;
						}
						else {
							IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
							if (containers.getAtomContainerCount()>0){
								containers = FileUtilities.parseSdf(inputFileName);					
								
//								int index = 0;
								for(IAtomContainer atc : containers.atomContainers()){
//									index ++;
									number_of_molecules++;
									System.out.println("\n\nMolecule no. " + number_of_molecules);
									try {
										 metabolites.add(mtf.findMetabolitesFromSequence(atc, 
													biotransformerSeqeuence, dmassesOrFormulas, massToleranceThreshold, annotate, opt));
										successful_predictions++;
									}
									catch(Exception e) {
										System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getLocalizedMessage());
									}
									
								}						
								
							}				
						}
					
						
						if(oFormat.contentEquals("csv")){
							FileUtilities.saveAtomContainerSetToCSV(metabolites, outputF);
						}
						else if(oFormat.contentEquals("sdf")){							
							SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputF));		
							sdfWriter.write(metabolites);
							sdfWriter.close();
						}					
						
					}
				}
				else{
					throw new IllegalArgumentException("\n\tFor metabolite identification, you must enter a list of masses, and/or molecular formulas");
				}
				
			}
			else if(task.contentEquals("pred")){
				if(biotransformerType != null) {
					if (biotransformerType.contentEquals("cyp450")){
						Cyp450BTransformer cyp450bt = new Cyp450BTransformer(BioSystemName.HUMAN);
						ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
						if (singleInput !=null){
							number_of_molecules++;
							biotransformations = cyp450bt.predictCyp450BiotransformationChain(singleInput, true, true, nrOfSteps, scoreThreshold);
							successful_predictions++;
						}
						else {
							IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
							if (containers.getAtomContainerCount()>0){
								containers = FileUtilities.parseSdf(inputFileName);					
								
//								int index = 0;
								for(IAtomContainer atc : containers.atomContainers()){
//									index ++;
									number_of_molecules++;
									System.out.println("\n\nMolecule no. " + number_of_molecules);
									try {
										biotransformations.addAll(cyp450bt.predictCyp450BiotransformationChain(atc, true, true, nrOfSteps, scoreThreshold));
										successful_predictions++;
									}
									catch(Exception e) {
										System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getLocalizedMessage());
									}
									
								}						
								
							}				
						}
										
						if(oFormat.contentEquals("csv")){
							cyp450bt.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
						}
						else if(oFormat.contentEquals("sdf")){
							cyp450bt.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
						}
					}
					
					else if (biotransformerType.contentEquals("ecbased")){
						ECBasedBTransformer ecbt =  new ECBasedBTransformer(BioSystemName.HUMAN);
			
						ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
						if (singleInput !=null){
							number_of_molecules++;
							biotransformations = ecbt.simulateECBasedMetabolismChain(singleInput, true, true, nrOfSteps, scoreThreshold);
							successful_predictions++;
						}
						else {
							IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
							if (containers.getAtomContainerCount()>0){
								containers = FileUtilities.parseSdf(inputFileName);
		//						biotransformations = ecbt.simulateECBasedMetabolismChain(containers, true, true, nrOfSteps, scoreThreshold);
								
//								int index = 0;
								for(IAtomContainer atc : containers.atomContainers()){
//									index++;
									number_of_molecules++;
									System.out.println("\n\nMolecule no. " + number_of_molecules);
		//							System.err.println("Molecule " + index );
									try {
										biotransformations.addAll(ecbt.simulateECBasedMetabolismChain(atc, true, true, nrOfSteps, scoreThreshold));
										successful_predictions++;
		//								throw new Exception();
									}
									catch(Exception e) {
										System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getMessage());
										continue;
									}
									
								}	
							
							}				
						}
										
						if(oFormat.contentEquals("csv")){
							ecbt.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
						}
						else if(oFormat.contentEquals("sdf")){
							ecbt.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
						}
					}
					else if (biotransformerType.contentEquals("hgut")){
						HGutBTransformer hgut = new HGutBTransformer();
		//				if(nrOfSteps!=8){
		//					System.out.println("\n=======>The number of steps for reductive metabolism is set. No need to set a number of steps for the human gut metabolism.\n\n");
		//				}
									
						ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
						if (singleInput !=null){
							number_of_molecules++;
							biotransformations = hgut.simulateGutMicrobialMetabolism(singleInput, true, true, nrOfSteps, scoreThreshold);
							successful_predictions++;
						}
						else {			
							IAtomContainerSet containers = FileUtilities.parseSdfAndAddTitles(inputFileName, hgut.inchiGenFactory);					
							if (containers.getAtomContainerCount()>0){
//								int index = 0;
								for(IAtomContainer atc : containers.atomContainers()){
//									index++;
									number_of_molecules++;
									System.out.println("\n\nMolecule no. " + number_of_molecules);
		//							System.err.println("Molecule " + index );
									try {
										biotransformations.addAll(hgut.simulateGutMicrobialMetabolism(atc, true, true, nrOfSteps, scoreThreshold));
										successful_predictions++;
		//								throw new Exception();
									}
									catch(Exception e) {
										System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getMessage());
										continue;
									}
									
								}	
							}				
						}
		//				System.out.println("No. of biotransformations: " + biotransformations.size());		
						if(oFormat.contentEquals("csv")){
							hgut.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
						}
						else if(oFormat.contentEquals("sdf")){
							hgut.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
						}						
					}
					else if (biotransformerType.contentEquals("phaseII")){
						Phase2BTransformer phase2b = new Phase2BTransformer(BioSystemName.HUMAN);
						
						ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
						if (singleInput !=null){
							number_of_molecules++;
							biotransformations = phase2b.applyPhase2TransformationsChainAndReturnBiotransformations(singleInput,
									true, true, true, nrOfSteps, scoreThreshold);
							successful_predictions++;
						}
						else {
							IAtomContainerSet containers = FileUtilities.parseSdfAndAddTitles(inputFileName, phase2b.inchiGenFactory);					
//							int index = 0;
							for(IAtomContainer atc : containers.atomContainers()){
//								index++;
								number_of_molecules++;
								System.out.println("\n\nMolecule no. " + number_of_molecules);
		//						System.err.println("Molecule " + index );
								try {
									biotransformations.addAll(phase2b.applyPhase2TransformationsChainAndReturnBiotransformations(atc, true, true, true, nrOfSteps, scoreThreshold));
									successful_predictions++;
		//							throw new Exception();
								}
								catch(Exception e) {
									System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getMessage());
									continue;
								}
								
							}
						
						}
										
						if(oFormat.contentEquals("csv")){
							phase2b.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
						}
						else if(oFormat.contentEquals("sdf")){
							phase2b.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
						}
					}
					
					else if (biotransformerType.contentEquals("superbio")){
						if(nrOfSteps!=12){
							System.out.println("\n\n=======>The configutration is set for this simulation. No need to set a number of steps for the super human transformer.\n\n");
						}
		
						HumanSuperBioTransformer hsbt = new HumanSuperBioTransformer();
						
						ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
						if (singleInput !=null){
							number_of_molecules++;
							biotransformations = hsbt.simulateHumanSuperbioMetabolism(singleInput,scoreThreshold);
							successful_predictions++;
						}
						else {
							IAtomContainerSet containers = FileUtilities.parseSdfAndAddTitles(inputFileName, hsbt.getInChIGenFactory());					
//							int index = 0;
							for(IAtomContainer atc : containers.atomContainers()){
//								index++;
								number_of_molecules++;
								System.out.println("\n\nMolecule no. " + number_of_molecules);
		//						System.err.println("Molecule " + index );
								try {
									biotransformations.addAll(hsbt.simulateHumanSuperbioMetabolism(atc,scoreThreshold));
									successful_predictions++;
		//							throw new Exception();
								}
								catch(Exception e) {
									System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getMessage());
									continue;
								}
								
							}
						
						}			
		
						if(oFormat.contentEquals("csv")){
							hsbt.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
						}
						else if(oFormat.contentEquals("sdf")){
							hsbt.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
						}
		
					}
					
					else if (biotransformerType.contentEquals("allHuman")){
						HumanSuperBioTransformer hsbt = new HumanSuperBioTransformer();
		
						ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
						if (singleInput !=null){
							number_of_molecules++;
							biotransformations = hsbt.predictAllHumanBiotransformationChain(singleInput, nrOfSteps, scoreThreshold);
							successful_predictions++;
						}
						else {
							IAtomContainerSet containers = FileUtilities.parseSdfAndAddTitles(inputFileName, hsbt.getInChIGenFactory());					
//							int index = 0;
							for(IAtomContainer atc : containers.atomContainers()){
//								index++;
								number_of_molecules++;
								System.out.println("\n\nMolecule no. " + number_of_molecules);
		//						System.err.println("Molecule " + index );
								try {
									biotransformations.addAll(hsbt.predictAllHumanBiotransformationChain(atc, nrOfSteps, scoreThreshold));
									successful_predictions++;
		//							throw new Exception();
								}
								catch(Exception e) {
									System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getMessage());
									continue;
								}
								
							}
						
						}			
		
						if(oFormat.contentEquals("csv")){
							hsbt.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
						}
						else if(oFormat.contentEquals("sdf")){
							hsbt.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
						}
						
		//				if (singleInput !=null){
		//					number_of_molecules++;
		//					if(oFormat.contentEquals("csv")){
		//						hsbt.predictAllHumanBiotransformationChainAndSaveToCSV(singleInput, nrOfSteps, scoreThreshold, outputF, annotate);
		//					}
		//					else if(oFormat.contentEquals("sdf")){
		//						hsbt.predictAllHumanBiotransformationChainAndSaveToSDF(singleInput, nrOfSteps, scoreThreshold, outputF, annotate);
		//					}
		//				}
		//				else {
		////					ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		//					IAtomContainerSet containers = FileUtilities.parseSdfAndAddTitles(inputFileName, hsbt.getInChIGenFactory());
		////					System.out.println(containers.getAtomContainerCount());
		//					number_of_molecules = containers.getAtomContainerCount();
		//					if(oFormat.contentEquals("csv")){
		//						hsbt.predictAllHumanBiotransformationChainAndSaveToCSV(containers, nrOfSteps, scoreThreshold, outputF, annotate);
		//					}
		//					else if(oFormat.contentEquals("sdf")){
		//						hsbt.predictAllHumanBiotransformationChainAndSaveToSDF(containers, nrOfSteps, scoreThreshold, outputF, annotate);
		//					}				
		//				}
					}
					
					else if (biotransformerType.contentEquals("env")){
						EnvMicroBTransformer ebt = new EnvMicroBTransformer();
						
						if (singleInput !=null){
							number_of_molecules++;
							if(oFormat.contentEquals("csv")){
		
								ebt.simulateEnvMicrobialDegradationAndSaveToCSV(singleInput, true, true, nrOfSteps, scoreThreshold, outputF, annotate);
								successful_predictions++;
							}
							else if(oFormat.contentEquals("sdf")){
								ebt.simulateEnvMicrobialDegradationAndSaveToSDF(singleInput, true, true, nrOfSteps, scoreThreshold, outputF, annotate);
								successful_predictions++;
							}
						}
						else {
							ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
							IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
							//biotransformations = ebt.applyEnvMicrobialTransformations(containers, true, true, scoreThreshold);
		//					biotransformations = ebt.simulateEnvMicrobialDegradation(containers, true, true, nrOfSteps, scoreThreshold);
							
//							int index = 0;
							for(IAtomContainer atc : containers.atomContainers()){
//								index++;
								number_of_molecules++;
								System.out.println("\n\nMolecule no. " + number_of_molecules);
		//						System.err.println("Molecule " + index );
								try {
									biotransformations.addAll(ebt.applyEnvMicrobialTransformationsChain(atc, true, true, nrOfSteps, scoreThreshold));
									successful_predictions++;
		//							throw new Exception();
								}
								catch(Exception e) {
									System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getMessage());
									continue;
								}
								
							}						
							
							if(oFormat.contentEquals("csv")){
								ebt.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
							}
							else if(oFormat.contentEquals("sdf")){
								ebt.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
							}					
						}
					
					
					}	
					
				}
				
				else if(biotransformerSeqeuence != null) {
					ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
					if (singleInput !=null){
						number_of_molecules++;
						biotransformations = biotransformerSeqeuence.runSequence(singleInput, scoreThreshold);
						successful_predictions++;
					}
					else {
						IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
						if (containers.getAtomContainerCount()>0){
							containers = FileUtilities.parseSdf(inputFileName);					
							
//							int index = 0;
							for(IAtomContainer atc : containers.atomContainers()){
//								index ++;
								number_of_molecules++;
								System.out.println("\n\nMolecule no. " + number_of_molecules + ": " + atc.getProperty(CDKConstants.TITLE));
								try {
									biotransformations.addAll(biotransformerSeqeuence.runSequence(atc, scoreThreshold));
									successful_predictions++;
								}
								catch(Exception e) {
									System.err.println("BioTransformer failed on molecule " + number_of_molecules + "\n" + e.getLocalizedMessage());
								}
								
							}						
							
						}				
					}
					UniversalBioTransformer hsbt = new UniversalBioTransformer();
					System.out.println("Format: " + oFormat);
					if(oFormat.contentEquals("csv")){
						System.out.println("saving to "+ outputF);
						hsbt.saveBioTransformationProductsToCSV(biotransformations, outputF, annotate);
					}
					else if(oFormat.contentEquals("sdf")){
						System.out.println("saving to "+ outputF);
						hsbt.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
					}
					
				}
				System.out.println("Successfully completed metabolism prediction for " + successful_predictions + " out of " + number_of_molecules + " molecule(s).");
			}
			else {
				throw new IllegalArgumentException("You entered an invalid biotransformer option.\n"
						+ "Choose one of the following: ecbased, cyp450, hgut, phaseII, superbio, allHuman, or env.\nType java -jar biotransformer-" + Version.current +".jar --help for help.");
			}			
		}		
	}

}




