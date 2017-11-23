package util;

import model.SequenceInfo;

import java.io.*;
import java.util.*;

/**
 * Created by nammi on 21/10/17.
 */
public class FileHelper {
    // HashMap()
    public static LinkedHashMap<String,String> readSequenceFile(String filename, LinkedHashMap<String,String> map){
        File file = new File(filename);
        FileReader fileReader = null;
        BufferedReader bufferedReader;
        try {
            fileReader = new FileReader(file);
        } catch (FileNotFoundException e){
            throw new IllegalArgumentException("File not found" + e);
        }
        bufferedReader = new BufferedReader(fileReader);
        String line;
        String id=null;
        StringBuffer stringBuffer = new StringBuffer();
        try {
            while ((line = bufferedReader.readLine()) != null) {
                if(line.contains(">")) {
                    if(id!=null){
                       // System.out.println(stringBuffer.toString());
                        map.put(stringBuffer.toString().toUpperCase(),id);
                        id =null;
                        stringBuffer.delete(0,stringBuffer.length());
                    }
                    String[] firstPart = line.split(">");
                    id = firstPart[1];
                    continue;
                }
                stringBuffer.append(line);
            }
            map.put(stringBuffer.toString().toUpperCase(),id);
        }catch (IOException e){
            throw new IllegalArgumentException(e.getMessage());
        }finally {
            try {
                if(bufferedReader!=null) {
                    bufferedReader.close();
                }
            }catch (IOException e){
                System.out.println("Warning Error while closing input stream");
            }
        }
        return map;
    }



    public static HashMap<Character, Integer> readAlphaFile(String fileName, HashMap<Character,Integer> letters){
        File file = new File(fileName);
        FileReader fileReader = null;
        BufferedReader bufferedReader;
        try {
            fileReader = new FileReader(file);
        } catch (FileNotFoundException e){
            throw new IllegalArgumentException("File not found" + e);
        }
        bufferedReader = new BufferedReader(fileReader);
        String line;
        StringBuffer stringBuffer = new StringBuffer();
        try {
            while ((line = bufferedReader.readLine()) != null) {
                stringBuffer.append(line);
            }
        }catch (IOException e){
            throw new IllegalArgumentException(e.getMessage());
        }
        finally {
            try {
                if(bufferedReader!=null) {
                    bufferedReader.close();
                }
            }catch (IOException e){
                System.out.println("Warning Error while closing input stream");
            }
        }
        char[] alpha = stringBuffer.toString().toUpperCase().toCharArray();
        for (int i=0;i<alpha.length;i++){
            letters.put(alpha[i],i);
        }
        return letters;
    }

    public static int[][] ReadScoringMatrix(String fileName, int[][] scoringMatrix){
        File file = new File(fileName);
        FileReader fileReader = null;
        BufferedReader bufferedReader;
        try {
            fileReader = new FileReader(file);
        } catch (FileNotFoundException e){
            throw new IllegalArgumentException("File not found" + e);
        }
        bufferedReader = new BufferedReader(fileReader);
        String line;
        StringBuffer stringBuffer = new StringBuffer();
        int i=0;
        int j=0;

        try {
            while ((line = bufferedReader.readLine()) != null) {
                StringTokenizer stringTokenizer = new StringTokenizer(line);
                while (stringTokenizer.hasMoreTokens()){
                    String token = stringTokenizer.nextToken();
                    try {
                        scoringMatrix[i][j] = Integer.parseInt(token.replaceAll("[^\\x00-\\x7F]", ""));
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                    j++;
                }
                i++;
                j=0;
            }
        }catch (IOException e){
            throw new IllegalArgumentException(e.getMessage());
        }
        finally {
            try {
                if(bufferedReader!=null) {
                    bufferedReader.close();
                }
            }catch (IOException e){
                System.out.println("Warning Error while closing input stream");
            }
        }
        return scoringMatrix;
    }
}
