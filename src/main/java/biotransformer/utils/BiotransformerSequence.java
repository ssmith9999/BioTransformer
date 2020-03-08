package biotransformer.utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.smiles.SmilesGenerator;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.btransformers.Cyp450BTransformer;
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.btransformers.EnvMicroBTransformer;
import biotransformer.btransformers.HGutBTransformer;
import biotransformer.btransformers.Phase2BTransformer;
import biotransformer.btransformers.Biotransformer.bType;
import biotransformer.transformation.Biotransformation;

public class BiotransformerSequence {
	//	public enum sequenceMode {
	//	
	//}
	
	protected LinkedHashMap<Biotransformer.bType, Integer> sequence;
	protected double scoreThreshold		= 0.0;
	protected Biotransformer utilityBiotransformer = null;
	

	public BiotransformerSequence(LinkedHashMap<Biotransformer.bType, Integer> mySequence, 
			double scoreThreshold) {
		this.sequence 		= mySequence;
		this.scoreThreshold = scoreThreshold;
	}	
	public BiotransformerSequence(LinkedHashMap<Biotransformer.bType, Integer> mySequence) {
		this.sequence = mySequence;
	}

	
	public BiotransformerSequence(String mySequence, double scoreThreshold, boolean annotate) throws Exception {
		this.sequence 		= createSequenceFromString(mySequence);
		this.scoreThreshold = scoreThreshold;
	}
	public BiotransformerSequence(String mySequence) throws Exception {
		this.sequence 		= createSequenceFromString(mySequence);
	}	
	
	private LinkedHashMap<Biotransformer.bType, Integer> createSequenceFromString(String mySequence) throws Exception{
		LinkedHashMap<Biotransformer.bType, Integer> seq = new LinkedHashMap<Biotransformer.bType, Integer>();
		String[] steps = mySequence.split("; ");
		for(String step : steps) {
			String[] attr = step.split(":");
//			System.out.println("mySequence: " + mySequence);
//			System.out.println("attr[0].toUpperCase(): " + attr[0].toUpperCase());
			seq.put(Biotransformer.bType.valueOf(attr[0].toUpperCase()), Integer.valueOf(attr[1]));
		}
		return seq;
	}
	
	public ArrayList<Biotransformation> runSequence(IAtomContainer startingCompound, double scoreThreshold) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		LinkedHashMap<Biotransformer.bType, Object> btransformers = new LinkedHashMap<Biotransformer.bType, Object>();
		IAtomContainerSet currentMetabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		currentMetabolites.addAtomContainer(startingCompound);
		
		for(Entry<bType, Integer> step : this.sequence.entrySet()) {
			
			ArrayList<Biotransformation> currentBiots =  new ArrayList<Biotransformation>();
			
			if(step.getKey() == bType.ALLHUMAN) {
				if(btransformers.get(step.getKey()) == null) {
					btransformers.put(step.getKey(), new HumanSuperBioTransformer());
				}
				currentBiots = ((HumanSuperBioTransformer) btransformers.get(step.getKey())).predictAllHumanBiotransformationChain(startingCompound, step.getValue(), this.scoreThreshold);
			}
			else if(step.getKey() == bType.CYP450) {
				if(btransformers.get(step.getKey()) == null) {
					btransformers.put(step.getKey(), new Cyp450BTransformer(BioSystemName.HUMAN));
				}
				currentBiots = ((Cyp450BTransformer) btransformers.get(step.getKey())).predictCyp450BiotransformationChain(startingCompound, true, true, step.getValue(), this.scoreThreshold);			
			}
			else if(step.getKey() == bType.ECBASED) {
				if(btransformers.get(step.getKey()) == null) {
					btransformers.put(step.getKey(), new ECBasedBTransformer(BioSystemName.HUMAN));
				}
				currentBiots = ((ECBasedBTransformer) btransformers.get(step.getKey())).simulateECBasedMetabolismChain(startingCompound, true, true, step.getValue(), this.scoreThreshold);		
			}
			else if(step.getKey() == bType.ENV) {
//				System.out.println("Starting ENV");
				if(btransformers.get(step.getKey()) == null) {
					btransformers.put(step.getKey(), new EnvMicroBTransformer());
				}
				currentBiots = ((EnvMicroBTransformer) btransformers.get(step.getKey())).applyEnvMicrobialTransformationsChain(startingCompound, true, true, step.getValue(), this.scoreThreshold);
//				System.out.println("Number of env. biotransformations: " + currentBiots.size());

			}			
			else if(step.getKey() == bType.HGUT) {
				if(btransformers.get(step.getKey()) == null) {
					btransformers.put(step.getKey(), new HGutBTransformer());
				}
				currentBiots = ((HGutBTransformer) btransformers.get(step.getKey())).simulateGutMicrobialMetabolism(startingCompound, true, true, step.getValue(), this.scoreThreshold);						
			}
			else if(step.getKey() == bType.PHASEII) {
				if(btransformers.get(step.getKey()) == null) {
					btransformers.put(step.getKey(), new Phase2BTransformer(BioSystemName.HUMAN));
				}
				currentBiots = ((Phase2BTransformer) btransformers.get(step.getKey())).applyPhase2TransformationsChainAndReturnBiotransformations(startingCompound, true, true, true, step.getValue(), this.scoreThreshold);										
			}
			else if(step.getKey() == bType.SUPERBIO) {
				throw new IllegalArgumentException("Invalid Argument: SUPERBIO cannot be used within a sequence, as it is already a customized sequence. Valid biotransformers within a sequence are ALLHUMAN, CYP450, ECBASED, ENVMICRO, HGUT, PHASEII.");
			}
			
			currentMetabolites.add(((Biotransformer) btransformers.get(step.getKey())).extractProductsFromBiotransformations(currentBiots));
			biotransformations.addAll(currentBiots);
		}
		
//		this.utilityBiotransformer = ((Collection<Entry<bType, Object>>) btransformers.entrySet()).stream().reduce((first, second) -> second).orElse(null)biotransformations.gtValue();
//		for(Biotransformation b: biotransformations) {
//			b.display();
//		}
		
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	

	public ArrayList<Biotransformation> runSequence(IAtomContainerSet startingCompounds, double scoreThreshold) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		SmilesGenerator smiGen 		= new SmilesGenerator().isomeric(); 
		for(IAtomContainer starting_ac : startingCompounds.atomContainers()) {
//			System.out.println(smiGen.create(starting_ac));
			try {
				biotransformations.addAll(runSequence(starting_ac, scoreThreshold));
			}
			catch (Exception e) {
				System.out.println(e);
			}
			
		}
		
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	
}
