
/*
 * CS4384.501
 * DFA Code
 * 
 * Team members:
 * Sebrina Kate Sendaydiego
 * Leonel Perez
 * Alexandria Andrade
 * Sirius Sharma
 * Chris Park
 */

import java.util.*;

public class DFA {

	private int currentState;
	private int[][] transitionTable;
	private int[] finalStates;
	private static String team = "Leonel Perez, Sebrina Kate Sendaydiego, Alexandria Andrade, Sirius Sharma, Chris Park";

	public DFA(int[][] transitionTable, int[] fOne) {
		this.currentState = 0;
		this.transitionTable = transitionTable;
		this.finalStates = fOne;
	}

	// setters and getters
	public int getCurrentState() {
		return currentState;
	}

	public void setCurrentState(int currentState) {
		this.currentState = currentState;
	}

	public int[] getFinalStates() {
		return finalStates;
	}

	public void setFinalStates(int[] finalStates) {
		this.finalStates = finalStates;
	}

	public int[][] getTransitionTable() {
		return transitionTable;
	}

	public void setTransitionTable(int[][] transitionTable) {
		this.transitionTable = transitionTable;
	}

	// returns the number of states
	public int getNumStates() {
		Set<Integer> states = new HashSet<>();
		int[][] tt = getTransitionTable();

		for (int i = 0; i < tt.length; i++) {
			states.add(i);

			for (int j = 0; j < tt[i].length; j++) {
				states.add(tt[i][j]);
			}
		}

		return states.size();
	}

	// parses through the string and uses transition table to get to the end state
	public boolean run(String input) {
		setCurrentState(0);

		for (int i = 0; i < input.length(); i++) {
			char inputSymbol = input.charAt(i);
			int inputIndex;
			if (inputSymbol == '0') {
				inputIndex = 0;
			} else {
				inputIndex = 1;
			}
			currentState = transitionTable[currentState][inputIndex];
		}

		// checks if end current state is one of the final states
		if (lookUp(currentState))
			return true;
		else
			return false;
	}

	public boolean lookUp(int state) {
		// System.out.println("End Current State: " + state);
		int[] finalsList = getFinalStates();

		// parses through the final states list to see if the state is valid
		for (int i = 0; i < finalsList.length; i++) {
			if (state == finalsList[i])
				return true;
		}
		return false;

	}

	// Returns Union of m1 and m2
	public static DFA union(DFA m1, DFA m2) {
		int[][] m1TT = m1.getTransitionTable();
		int[][] m2TT = m2.getTransitionTable();
		int statesNum = m1.getTransitionTable().length * m2.getTransitionTable().length;
		int[][] tInt = new int[statesNum][2]; // transition table for union DFA

		List<Integer> finalsList = new ArrayList<Integer>(); // list the final states of the union DFA

		// fill transition table for union DFA
		for (int i = 0; i < statesNum; i++) {
			int state1 = i / m2TT.length;
			int state2 = i % m2TT.length;

			for (int j = 0; j < 2; j++) {
				int next1 = m1TT[state1][j];
				int next2 = m2TT[state2][j];
				int next = next1 * m2TT.length + next2;
				tInt[i][j] = next;
			}

			// determine if the state is final and add to final states list
			if (m1.lookUp(state1) || m2.lookUp(state2)) {
				finalsList.add(i);
			}
		}

		// create array of all final states
		int[] finalStates = new int[finalsList.size()];
		for (int i = 0; i < finalsList.size(); i++) {
			finalStates[i] = finalsList.get(i);
		}

		// intersection of the two DFAs
		DFA unionDFA = new DFA(tInt, finalStates);

		return unionDFA;
	}

	// Returns Intersection of m1 and m2
	public static DFA intersection(DFA m1, DFA m2) {

		int[][] m1TT = m1.getTransitionTable();
		int[][] m2TT = m2.getTransitionTable();
		int statesNum = m1.getTransitionTable().length * m2.getTransitionTable().length;
		int[][] tInt = new int[statesNum][2]; // transition table for union DFA

		List<Integer> finalsList = new ArrayList<Integer>(); // list the final states of the union DFA

		// fill transition table for union DFA
		for (int i = 0; i < statesNum; i++) {
			int state1 = i / m2TT.length;
			int state2 = i % m2TT.length;

			for (int j = 0; j < 2; j++) {
				int next1 = m1TT[state1][j];
				int next2 = m2TT[state2][j];
				int next = next1 * m2TT.length + next2;
				tInt[i][j] = next;
			}

			// determine if the state is final and add to final states list
			if (m1.lookUp(state1) && m2.lookUp(state2)) {
				finalsList.add(i);
			}
		}

		// create array of all final states
		int[] finalStates = new int[finalsList.size()];
		for (int i = 0; i < finalsList.size(); i++) {
			finalStates[i] = finalsList.get(i);
		}

		// intersection of the two DFAs
		DFA intersectDFA = new DFA(tInt, finalStates);

		return intersectDFA;
	}

//	// Returns Difference of m1 and m2 (order matters: L(m1)-L(m2))
	public static DFA difference(DFA m1, DFA m2) {
		DFA difDFA = DFA.intersection(m1, DFA.complement(m2));

		return difDFA;
	}

	// Returns Complement of m
	public static DFA complement(DFA m) {
		// transition table stays the same
		// final states are swapped w/ non final

		int maxStates = m.getTransitionTable().length;
		int[] oldFinArr = m.getFinalStates();

		List<Integer> oldFinList = new ArrayList<Integer>(); // contains old final states for m
		List<Integer> newFinList = new ArrayList<Integer>(); // contains new final states for complement

		// fill list to use .contains later
		for (int i = 0; i < oldFinArr.length; i++) {
			oldFinList.add(oldFinArr[i]);
		}

		// if cur state is not a final state, add it to new list
		for (int state = 0; state < maxStates; state++) {
			if (!oldFinList.contains(state))
				newFinList.add(state);
		}

		// convert new list to a final state array
		int[] finalStates = new int[newFinList.size()];

		for (int i = 0; i < newFinList.size(); i++) {
			finalStates[i] = newFinList.get(i);
		}

		// new DFA complement of m
		DFA compDFA = new DFA(m.getTransitionTable(), finalStates);

		return compDFA;
	}

	// Returns String with names of all group members
	public static String credits() {
		return team;
	}

	// Test if the only lnaguage accepted bya dfa is Empty set
	public boolean isEmptyLanguage() {
		Set<Integer> finalStateSet = new HashSet<>();
		Set<Integer> visited = new HashSet<>();
		Queue<Integer> queue = new LinkedList<>();
		queue.add(0); // initial state

		for (int state : finalStates) {
			finalStateSet.add(state);
		}

		while (!queue.isEmpty()) {
			currentState = queue.poll();
			visited.add(currentState);

			// if current state is final, DFA does not accept only empty
			if (finalStateSet.contains(currentState))
				return false;

			// add all next states to queue
			for (int i = 0; i < transitionTable[currentState].length; i++) {
				int nextState = transitionTable[currentState][i];

				if (!visited.contains(nextState)) {
					queue.add(nextState);
				}
			}
		}

		// if we have visited all states and none of them are final
		return true;

	}

	// Test if the language accepted by a DFA is every string
	public boolean isUniversalLanguage() {
		Set<Integer> finalStateSet = new HashSet<>();
		Set<Integer> visited = new HashSet<>();
		Queue<Integer> queue = new LinkedList<>();
		int finalCounter = 0;
		int stateCounter = 0;
		queue.add(0); // initial state

		for (int state : finalStates) {
			stateCounter++;
			finalStateSet.add(state);
		}

		while (!queue.isEmpty()) {
			currentState = queue.poll();
			visited.add(currentState);

			// if current state is final, DFA does not accept only empty
			if (finalStateSet.contains(currentState))
				finalCounter++;

			// add all next states to queue
			for (int i = 0; i < transitionTable[currentState].length; i++) {
				stateCounter++;
				int nextState = transitionTable[currentState][i];

				if (!visited.contains(nextState)) {
					stateCounter++;
					queue.add(nextState);
				}
			}
		}

		// if we have visited all states and none of them are final
		// System.out.println("States vs. Final states: " + stateCounter + " / " +
		// finalCounter);
		if (stateCounter == finalCounter) {
			return true;
		} else {
			return false;
		}

	}

	// checks if the DFA is a subset of m2
	public boolean isSubsetOf(DFA m2) {
		DFA thisDFA = new DFA(getTransitionTable(), getFinalStates());

		if (!difference(thisDFA, m2).isEmptyLanguage())
			return false;

		return true;
	}

	// checks if language is infinite
	public boolean isInfinite() {
		int length = getNumStates();
		String bitString = "";

		// loops thru every possible string in the range of length n to 2n
		while (length < 2 * getNumStates()) {
			for (int i = 0; i < Math.pow(2, length); i++) {
				bitString = createBit(i, length);

				if (run(bitString))
					return true;
			}
			length++;
		}
		return false;
	}

	// converts value into binary and length of the bit string
	public static String createBit(int value, int length) {
		String bitStr = "";
		BitSet bitSet = BitSet.valueOf(new long[] { value });

		for (int i = length - 1; i >= 0; i--) {
			bitStr += bitSet.get(i) ? "1" : "0";
		}

		return bitStr;
	}

	// Checks if L(m1) = L(m2)
	public boolean equals(DFA m2) {
		DFA m1 = new DFA(getTransitionTable(), getFinalStates());

		// if L(m1)-L(m2) and L(m2)-L(m1) is null, then L(m1) = L(m2)
		if (difference(m1, m2).isEmptyLanguage() && difference(m2, m1).isEmptyLanguage())
			return true;

		return false;
	}

	// Converts DFA to bit string
	public static String compress(DFA dfa1) {
		StringBuilder bitString = new StringBuilder();
		int[][] tt1 = dfa1.getTransitionTable();
		int[] finalStates = dfa1.getFinalStates();
		int numStates = tt1.length;
		int numInputs = tt1[0].length;
		int numFinalStates = finalStates.length;

		// encode the number of states, inputs, and final states into the string
		bitString.append(numStates);
		bitString.append(numInputs);
		bitString.append(numFinalStates);

		// encode the transition table
		for (int[] row : tt1) { // goes through each transition table row
			for (int num : row) { // converts the integer values at each index into binary string
				String binaryString = Integer.toBinaryString(num);
				while (binaryString.length() < 3) {
					binaryString = "0" + binaryString;
				}
				bitString.append(binaryString); // appends binary strings together
			}
		}

		// encode the final states
		for (int finalState : finalStates) {
			String binaryString = Integer.toBinaryString(finalState);
			while (binaryString.length() < 2) {
				binaryString = "0" + binaryString;
			}
			bitString.append(binaryString); // append binary string of final states
		}

		// System.out.println("Bit String from given DFA: " + bitString.toString());
		return bitString.toString();
	}

	// Converts bit string to DFA
	public static DFA decompress(String compressed) {
		int numStates = compressed.charAt(0) - '0'; // get number of states
		int numInputs = compressed.charAt(1) - '0'; // get number of inputs
		int numFinalStates = compressed.charAt(2) - '0'; // get number of final states
		int tableBegins = 3; // starting index of transition table
		int[] finalStates = new int[numFinalStates]; // stores final states

		// get final states from compressed string
		for (int i = 0; i < numFinalStates; i++) {
			finalStates[i] = compressed.charAt(3 + i) - '0';
		}

		int[][] transitions = new int[numStates][numInputs]; // stores transition table

		// get transition table from compressed string
		for (int i = 0; i < numStates; i++) {
			for (int j = 0; j < numInputs; j++) {
				int startIndex = tableBegins + (i * numInputs + j) * 3; // find starting index for transition table
																		// representation
				String binaryString = compressed.substring(startIndex, startIndex + 3);
				int transition = Integer.parseInt(binaryString, 2);
				transitions[i][j] = transition;
			}
		}
		// construct and return DFA
		return new DFA(transitions, finalStates);
	}

	// checks if the new DFA is identical to the original DFA
	public boolean identical(DFA dfa2) {
		String compressed1 = compress(this);
		String compressed2 = compress(dfa2);
		DFA decompressed1 = decompress(compressed1);
		DFA decompressed2 = decompress(compressed2);

		int[][] table1 = decompressed1.getTransitionTable();
		int[][] table2 = decompressed2.getTransitionTable();
		int[] finalStates1 = decompressed1.getFinalStates();
		int[] finalStates2 = decompressed2.getFinalStates();

		// check if the final state or transition table lengths are not equal
		if (finalStates1.length != finalStates2.length || table1.length != table2.length) {
			return false;
		}

		// compare transition tables
		for (int i = 0; i < table1.length; i++) {
			for (int j = 0; j < table1[i].length; j++) {
				int value1 = table1[i][j];
				int value2 = table2[i][j];

				if (value1 != value2) {
					return false;
				}
			}
		}

		// compare final states
		for (int i : finalStates2) {
			if (!decompressed1.lookUp(i)) {
				return false;
			}
		}

		return true;
	}
}
