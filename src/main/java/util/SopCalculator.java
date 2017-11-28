package util;

import model.Sequence;
import model.SequenceInfo;

import java.util.*;

/**
 * Created by nammi on 28/11/17.
 */
public class SopCalculator {

    public static int[][] populateScoringMatrix(int[] scoringArray,int[][] scoringmatrix,char[] letters)
    {
        int i,j;
        int m_pos = 0;
        m_pos = 0;
        for (i = 0;i < letters.length;i++){
            for (j = 0;j <= i;j++){
                if (i == j){
                    scoringmatrix[letters[i]-65][letters[j]-65] += scoringArray[m_pos];
                }else{
                    scoringmatrix[letters[i]-65][letters[j]-65] += scoringArray[m_pos];
                    scoringmatrix[letters[j]-65][letters[i]-65] += scoringArray[m_pos];
                }
                m_pos++;
            }
        }
        return scoringmatrix;
    }

   public static int calcSop(String filename, int[][] scoringMatrix,int gap){
        LinkedHashMap<String,String> sequences = new LinkedHashMap<>();
        sequences = FileHelper.readSequenceFile(filename, sequences);
        String[] seq = sequences.keySet().toArray(new String[sequences.size()]);

        int totalScore =0;
        int count =0;
        for (int i=0;i<seq.length;i++)
        {
            for(int j= i+1;j<seq.length;j++){
               totalScore+= calulateSP(seq[i],seq[j],scoringMatrix,gap);
                count++;
            }
        }

        return totalScore/count;

    }

    public static int calulateSP(String seq1, String seq2,int[][] scroringMatrix, int gap){

        int score=0;
        for(int i=0;i<seq1.length();i++){
            if(!((seq1.charAt(i)=='.')&&(seq2.charAt(i) == '.'))){
                if(seq1.charAt(i)=='.'||seq2.charAt(i)=='.'){
                    score+=gap;
                }else {
                    System.out.println(seq1.charAt(i) + " " +seq2.charAt(i) );
                    score+=scroringMatrix[seq1.charAt(i)-65][seq2.charAt(i)-65];
                }
            }
        }

    return score;
    }
}
