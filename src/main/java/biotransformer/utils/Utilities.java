package biotransformer.utils;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Map.Entry;

import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import biotransformer.transformation.Biotransformation;

public class Utilities {

	public Utilities() {
		// TODO Auto-generated constructor stub
	
	}
	
	public static ArrayList<String> removeDuplicateStrings(ArrayList<String> listOfStrings){
		ArrayList<String> unique =  new ArrayList<String>();		
		LinkedHashSet<String> set = new LinkedHashSet<String>(listOfStrings);
		unique = new ArrayList<String>(set);
		
		return unique;
	}
	
	public static void print(ArrayList<String> aList){
		for(String i : aList){
			System.out.println(i);
		}
	}

	public static IAtomContainerSet getCDKAtomContainerSet(){
		return DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		
	}
	
	public static void annotateAtomContainerWithProps(IAtomContainer molecule, LinkedHashMap<Object, Object> props){
		LinkedHashMap<Object, Object> p = new LinkedHashMap<Object, Object>();
		for(Entry<Object, Object> m: props.entrySet()){
			if(m.getKey() != "Synonym"){
				p.put(m.getKey(), m.getValue());
			}
			else {
				p.put(CDKConstants.TITLE,m.getValue());
			}
		}		
//		System.out.println(p);
		molecule.setProperties(p);	
	}


	public static ArrayList<Biotransformation> selectUniqueBiotransformations(ArrayList<Biotransformation> biotransformations){
		ArrayList<Biotransformation> unique_bts = new ArrayList<Biotransformation>();
		for(int i = 0; i < biotransformations.size(); i ++){
			if(!containsBiotransformation(unique_bts, biotransformations.get(i))){
				unique_bts.add(biotransformations.get(i));
			}
		}	
		return unique_bts;
	}
	
	public static boolean containsBiotransformation(ArrayList<Biotransformation> biotransformations, Biotransformation bt){
		boolean inc = false;
		for(int i = 0; i < biotransformations.size(); i ++){
			if(biotransformations.get(i).equals(bt)){
				inc = true;
				break;
			}
		}
		
		
		return inc;
	}
	public static IAtomContainerSet createEmptyAtomContainerSet(){
		return DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
	}
	
//	public static IAtomContainerSet uniquefy(IAtomContainerSet molecules)
//			throws Exception {
//		if (molecules != null && (!molecules.isEmpty()) && molecules.getAtomContainerCount() > 1) {
//			
//			IAtomContainerSet uniqueContainer = DefaultChemObjectBuilder.getInstance()
//					.newInstance(IAtomContainerSet.class);
//			
//			uniqueContainer.addAtomContainer(molecules.getAtomContainer(0));
//			
//			for (int i = 1; i < molecules.getAtomContainerCount(); i++) {
//				if (! ( (molecules.getAtomContainer(i) == null) || atomContainerInclusionHolds(uniqueContainer,
//						molecules.getAtomContainer(i) ))) {
//					uniqueContainer.addAtomContainer(molecules.getAtomContainer(i));
//				}
//			}
//
//			return uniqueContainer;
//		}
//
//		else
//			return molecules;
//
//	}

}
