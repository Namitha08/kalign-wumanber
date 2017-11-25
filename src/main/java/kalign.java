import model.SequenceInfo;
import util.FileHelper;
import util.MatrixHelper;

import java.util.*;
/*Comment 1*/

/**
 * Created by nammi on 22/11/17.
 */
public class kalign {

    public static int numseq;
    public static int numprofiles;

    public static void main(String[] args) {
        int[][][] matches = new int[8000][][];
        double[][] dm;
        int a,b;
        int[] tree;

        SequenceInfo si = new SequenceInfo();
        si = readSequence("input.fasta", si);
        for (int i = 8000-1;i>=0;i--){
            matches[i] = null;
        }
        fill_hash(si,matches);
        dm = new double[numprofiles][];
        for (int i = numprofiles-1;i>=0;i--){
            dm[i] = new double[numprofiles];
            for (int j = numprofiles-1;j>=0;j--){
                dm[i][j] = 0;
            }
        }
        System.out.println("Distance Calculation:");
        b = (numseq*(numseq-1))/2;
        a = 1;
        for (int i = 0; i < numseq;i++){
                for (int j = i+1; j < numseq;j++){
                    dm[i][j] = distance_calculation(matches,si.sl[i],si.sl[j],i,j);
                    System.out.println("percent done" + (double)a /(double)b * 100);
                    a++;
                }
            }

            MatrixHelper.printMatrix(dm);

        for (int i = 8000-1; i>=0;i--){
            if (matches[i] !=null){
                for (int j = numseq-1;j>=0;j--){
                    if (matches[i][j]!=null){
                        matches[i][j][0] &= 0x0000ffff;
                    }
                }
            }
        }
        tree = new int[(numseq-1)*2];
        tree = upgma(dm,tree);
        MatrixHelper.printArray(tree);
        MatrixHelper.printMatrix(dm);
    }

    static int[] upgma(double[][]dm,int[] tree)
    {
        int i,j,t;
	    int[] as = new int[numprofiles];
        double max;
        int node_a = 0;
        int node_b = 0;
        int cnode;
        cnode = numseq;

        for (i = numseq-1;i>=0; i--){
            as[i] = 1;
        }
        for (i = numseq; i < numprofiles;i++){
            as[i] = 0;
        }
        t = 0;
        while (cnode != numprofiles){
            max = Integer.MIN_VALUE;
            for (i = 0;i < numprofiles; i++){
                if (as[i]!=0){
                    for ( j = i +1;j < numprofiles;j++){
                        if (as[j]!=0){
                            if (dm[i][j] > max){
                                max = dm[i][j];
                                node_a = i;
                                node_b = j;
                            }
                        }
                    }
                }
            }
		/*deactivate  sequences to be joined*/
            as[node_a] = 0;
            as[node_b] = 0;
            tree[t] = node_a;
            tree[t+1] = node_b;
		/*calculate new distances*/
            for (i = numprofiles-1;i>=0;i--){
                if (as[i]!=0){
                    dm[i][cnode] = (node_a < i)?dm[node_a][i]:dm[i][node_a];
                    dm[i][cnode] += (node_b < i)?dm[node_b][i]:dm[i][node_b];
                    dm[i][cnode] /= 2;
                }
            }
            as[cnode] = 1;
            cnode++;
            t += 2;
        }
        return tree;
    }

    static double distance_calculation(int[][][] matches,int len_a,int len_b,int a,int b)//,int size)
    {
        double d = 0;
        int i,j,c,tmp1,tmp2;
	    int[] diagonal;
	    int[] p1;
	    int[] p2;
	    int[] p=null;
        diagonal = new int[len_a+len_b];
        for(i = len_a + len_b-1;i>=0;i--){
            diagonal[i] = 0;
        }
      // int temp = diagonal[len_a];
        for (c = 8000-1;c>=0;c--){
            if (matches[c] !=null){
                if (((matches[c][a][0]&0x0000ffff)-1) !=0){
                    p1 = matches[c][a];
                    tmp1 = (p1[0] >> 16);
                    if (((matches[c][b][0]&0x0000ffff)-1)!=0){
                        p2 = matches[c][b];
                        tmp2 = tmp1 * (p2[0]>>16);
                        System.out.println("Printing p's");
                        for (i = (p1[0] & 0x0000ffff)-1;i>0;i--){
                           // int temp_p = p[temp+p1[i]];
                           // p = diagonal - p1[i];
                            for (j = (p2[0] & 0x0000ffff)-1;j>0;j--){
                                //TODO :  check what shld be p ?
                                diagonal[len_a-p1[i] + p2[j]] += tmp2;
                                System.out.print( diagonal[len_a-p1[i] + p2[j]] + " ");
                            }
                            System.out.println("\n");
                        }
                    }
                }
            }
        }

      //  diagonal -= len_a;
        c = 0; //min
        for (i = 3-1; i>=0;i--){
            for(j = len_a+len_b-1;j>=0;j--){
                //fprintf(stdout,"	%d	c:%d\n",diagonal[j],c);
                if (diagonal[j]!=0){
                    if(diagonal[j] > (c & 0x0000ffff)){
                        c = diagonal[j]| (j<<16);
                        //fprintf(stderr,"c:%d	%d\n",c &0x0000ffff, c >>16);
                    }
                }
            }
            if (c!=0){
                d += c & 0x0000ffff;
                diagonal[c>>16] = 0;
                c = 0;
            }
        }
       // free(diagonal);
        return d;
    }
    public static SequenceInfo readSequence(String queryFile, SequenceInfo si){

        LinkedHashMap<String,String> sequences = new LinkedHashMap<String, String>();
        sequences = FileHelper.readSequenceFile(queryFile, sequences);
        Iterator iterator =  sequences.entrySet().iterator();
        int i=0;
        numseq = sequences.size();
        numprofiles = (numseq *2 ) - 1;
        si.s =  new int[numseq][];
        si.sn =  new String[sequences.size()];
        si.sl = new int[numprofiles];
        si.lsn =  new int[numseq];
        while (iterator.hasNext()){
            Map.Entry pair = (Map.Entry) iterator.next();
            //TODO : Extra last space not allocated
            String key = (String) pair.getKey();
            String value = (String) pair.getValue();
            String seq =  key.toUpperCase();
            // TODO : Allocating 1 more letter space
            si.s[i] = new int[key.length()+1];
            for (int j=0;j<key.length();j++){
                si.s[i][j] = seq.charAt(j)-65;
            }
            // TODO : assigning 0 to last extra space
            si.s[i][key.length()] = 0;
            MatrixHelper.printArray(si.s[i]);
            si.sl[i] = key.length();
            si.sn[i] = value;
            si.lsn[i] = value.length();
            i++;
        }

        int c = 0;
        int n = 0;
        i = 0;
        int j = 0;

        si.sip = new int[numprofiles][];
        si.nsip = new int[numprofiles];
        si.gis = new int[numprofiles][];
        si.relpos = new int[numprofiles][];

        for(i= numseq-1 ;i>=0; i--) {
            //TODO : Extra last space allocated
            si.gis[i] = new int[si.sl[i] + 1];
            si.relpos[i] = new int[si.sl[i] + 1];


            for (j = si.sl[i]; j >= 0; j--) {
                si.gis[i][j] = 0;
                si.relpos[i][j] = j;
            }
        }
        return si;
    }

    public static void fill_hash(SequenceInfo si,int matches[][][])
    {
        int i,j,c,f;
        int key = 0;
        // size 26
	    int aacode[] = {0,0,1,2,3,4,5,6,7,-1,8,9,10,11,-1,12,13,14,15,16,-1,17,18,0,19,0};
        // size 20
        int aadecode[] = {0,2,3,4,5,6,7,8,10,11,12,13,15,16,17,18,19,21,22,24};
        //size 30
        int patterns[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        // size 10
	    int map[] = {0,0,0,0,0,0,0,0,0,0};
        int[] sp;
        int out[][] = null;
        System.out.println("Printing Keys\n");
        for (i = numseq-1;i>=0;i--){
            sp = si.s[i];
            for(j = si.sl[i]-2-1;j>=0;j--){
                key = ((aacode[sp[j]]*20)+aacode[sp[j+1]])*20+aacode[sp[j+2]];
                //System.out.print(key+" ");
                if (matches[key]==null){
                   // System.out.println("key "+ key);
                    matches[key] = new int [numprofiles][];
                    for (c = numprofiles-1;c>=0;c--){
                        matches[key][c] = null;
                    }
                }
                if (matches[key][i] == null){
                    System.out.println("key " + key + " i " + i + " j " + j);
                    matches[key][i] = new int[10];
                    matches[key][i][0] = 1;
                }
                if ((matches[key][i][0] %10) == 0){
                   System.out.println("Realloc key " + key + " i " + i + " j " + j);
                    matches[key][i] = Arrays.copyOf(matches[key][i], (matches[key][i][0] + 10));
                }
                matches[key][i][matches[key][i][0]] = j;
                matches[key][i][0] += 1;
            }
            System.out.println("\n");
        }

        for (i = 8000-1;i>=0;i--){
            if (matches[i] !=null){
                for (j = numseq-1;j>=0;j--){
                    if (matches[i][j] != null){
                        matches[i][j][0] |= 0x00040000;
                    }
                }
            }
        }

        c = 0;
        for (i = numseq-1;i>=0;i--){
            for (j = 8000-1;j>=0;j--){
                if (matches[j] != null){
                   // System.out.println("out j " + j);
                    if(matches[j][i] == null){
                       // System.out.println("j " + j + " i " + i +"\n");
                        map[c] = j;
                        key = j;
                        patterns[c*3 + 2] = aadecode[key / 400];
                        key %= 400;
                        patterns[c*3 + 1] = aadecode[key /20];
                        key %= 20;
                        patterns[c*3 + 0] = aadecode[key];
                        c++;
                    }
                    if(c == 10){
                        //start of 10x Wu-Manber;
                        out = ten_wu_manber(si.s[i],si.sl[i],patterns);
                        System.out.println("c = 10 ");
                        MatrixHelper.printMatrix(out);
                        for (f = 0;f < 10;f++){
                            matches[map[f]][i] = out[f];
                            matches[map[f]][i][0] |= 0x00010000;
                        }
                        //free(out);
                        c = 0;
                    }
                }
            }
            if (c!=0){
                for (f = c*3;f < 30;f++){
                    patterns[f] = 9;
                }
                out = ten_wu_manber(si.s[i],si.sl[i],patterns);
                System.out.println("c!=0");
                MatrixHelper.printMatrix(out);
                for (f = 0;f < c;f++){
                    matches[map[f]][i] = out[f];
                    matches[map[f]][i][0] |= 0x00010000;
                }
                //free(out);
                c = 0;
            }
        }
        MatrixHelper.printArray(patterns);

    }

    public static int[][] ten_wu_manber(int[] seq,int len,int p[])
    {
        // unsigned
        int Tc = 0;
        int Tc2 = 0;
        // register
        int R0 = 0;
        int R1= 153391689;//1;
        int S0 = 0;
        int S1 = 0;
        // size 26
        int T[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        // size 26
        int T2[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	    int [][]out = null;
        out = new int[10][];
        for (Tc = 10-1;Tc>=0;Tc--){
            out[Tc] = new int[10];
            out[Tc][0] = 1;
        }

        for (Tc = 0; Tc < 30;Tc++){
            T[p[Tc]] |= 1 << Tc;
            T2[p[Tc]] = T[p[Tc]]&153391689;//MASK;
        }
        while(len>0){
            //fprintf(stderr,"%d\n",len);
            len--;
            Tc = T[seq[len]];
            Tc2 = T2[seq[len]];
            S0 <<= 1;
            S0 |= 153391689;//MASK;
            S0 &= Tc;
            S0 |= Tc2;
            //S1 = ((R1 << 1) & Tc) | ((R0 | S0) << 1) | R0;//| 1;
            S1 = (R1 << 1) & Tc;//| 1;
            S1 |= (R0 | S0) << 1;
            S1 |= R0;
            S1 |= 153391689;

            //S1 = (R1 << 1) & (Tc) | ((R0 | S0) << 1) | (R0) | (153391689);

            if ((S1 & 0x24924924) !=0 ){
                if ((S1 & 0x00004924) != 0 ){
                    if ((S1 & 0x00000024) !=0){
                        if ((S1 & 4) !=0){
                            if ((out[0][0] %10) !=0){
                                out[0] = Arrays.copyOf(out[0],(out[0][0]+10));
                            }
                            out[0][out[0][0]] = len;
                            out[0][0]++;
                        }
                        if ((S1 & 32) !=0){
                            if ((out[1][0] %10) !=0){
                                out[1] = Arrays.copyOf(out[1],(out[1][0]+10));
                            }
                            out[1][out[1][0]] = len;
                            out[1][0]++;
                        }
                    }
                    if ((S1 & 0x00004900) !=0){
                        if ((S1 & 256)!=0){
                            if ((out[2][0] %10) !=0){
                                out[2] = Arrays.copyOf(out[2],(out[2][0]+10));
                            }
                            out[2][out[2][0]] = len;
                            out[2][0]++;
                        }
                        if ((S1 & 2048) !=0){
                            if ((out[3][0] %10) !=0){
                                out[3] = Arrays.copyOf(out[3],(out[3][0]+10));
                            }
                            out[3][out[3][0]] = len;
                            out[3][0]++;
                        }
                        if ((S1 & 16384)!=0){
                            if ((out[4][0] %10) !=0){
                                out[4] = Arrays.copyOf(out[4],(out[4][0]+10));
                            }
                            out[4][out[4][0]] = len;
                            out[4][0]++;
                        }
                    }
                }



                if ((S1 & 0x24920000) !=0){
                    if ((S1 & 0x00920000) !=0){
                        if ((S1 & 131072) !=0){
                            if ((out[5][0] %10) !=0){
                                out[5] = Arrays.copyOf(out[5],(out[5][0]+10));
                            }
                            out[5][out[5][0]] = len;
                            out[5][0]++;
                        }
                        if ((S1 & 1048576) !=0){
                            if ((out[6][0] %10) !=0){
                                out[6] = Arrays.copyOf(out[6],(out[6][0]+10));
                            }
                            out[6][out[6][0]] = len;
                            out[6][0]++;
                        }
                        if ((S1 & 8388608) !=0){
                            if ((out[7][0] %10) !=0){
                                out[7] = Arrays.copyOf(out[7],(out[7][0]+10));
                            }
                            out[7][out[7][0]] = len;
                            out[7][0]++;
                        }
                    }

                    if ((S1 & 0x24000000) !=0){
                        if ((S1 & 67108864) !=0){
                            if ((out[8][0] %10) !=0){
                                out[8] = Arrays.copyOf(out[8],(out[8][0]+10));
                            }
                            out[8][out[8][0]] = len;
                            out[8][0]++;
                        }
                        if ((S1 & 536870912) !=0){
                            if ((out[9][0] %10) !=0){
                                out[9] = Arrays.copyOf(out[9],(out[9][0]+10));
                            }
                            out[9][out[9][0]] = len;
                            out[9][0]++;
                        }
                    }
                }
            }
            R0 = S0;
            R1 = S1;
        }
        return out;
    }
}


