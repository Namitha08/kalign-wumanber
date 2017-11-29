import model.DPStructure;
import model.SequenceInfo;
import util.FileHelper;
import util.MatrixHelper;
import util.SopCalculator;

import java.util.*;
/*Comment 1*/

/**
 * Created by nammi on 22/11/17.
 */
public class kalign {

    public static int numseq;
    public static int numprofiles;
	public static int INFINITY = 0x8000000;
	static int gpo = 61;
	static int gpe = 18;

//	public static void main(String[] args) {
//
//		int[][] scoringMatrix = new int[26][26];
//		int[] scoringArray = new int[276];
//		char[] letters = {'A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'};
//		scoringArray = FileHelper.ReadScoringArray("gon250.txt",scoringArray);
//		for(int i=0;i<scoringMatrix.length;i++){
//			for (int j=0;j<scoringMatrix[i].length;j++){
//				scoringMatrix[i][j]=0;
//			}
//		}
//
//		scoringMatrix = SopCalculator.populateScoringMatrix(scoringArray,scoringMatrix,letters);
//		int score = SopCalculator.calcSop("out.fasta",scoringMatrix,-61);
//		System.out.println(score);
//	}
//
    public static void main(String[] args) {
        int dia = 24;
        int[][][] matches = new int[8000][][];
        double[][] dm;
        int a,b;
        int[] tree;
        int scoringArray[] =  new int[276];
        int[][] submatrix = new int[32][32];
        int len_a;
        int len_b;
        int path[];
        int profa[];
        int profb[];
        char[] letters = {'A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'};
        SequenceInfo si = new SequenceInfo();
        si = readSequence("input.fasta", si);
        DPStructure dp = new DPStructure();
//        for (int i = 8000-1;i>=0;i--){
//            matches[i] = null;
//        }
        int[][] profile = new int[numprofiles][];
		int newnode = numseq;
//        for (int i = numprofiles-1; i>=0;i--){
//            profile[i] = null;
//        }
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

           // MatrixHelper.printMatrix(dm);

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
//        MatrixHelper.printMatrix(dm);

        scoringArray = FileHelper.ReadScoringArray("gon250.txt", scoringArray);

    //   MatrixHelper.printArray(scoringArray);

        for (int i = 32-1;i>=0 ;i--){
            for (int j = 32-1;j>=0;j--){
                submatrix[i][j] = gpe;
            }
        }

        submatrix = populateScoringMatrix(scoringArray, submatrix, gpo, letters);
        //MatrixHelper.printMatrix(submatrix);

        dp = initializeDPStructure(dp,511,511);

        int j;
		MatrixHelper.printArray(si.sl);
        System.out.println("\nAlignment:\n");
        for (int i = 0; i < (numseq-1)*2;i +=2){
            System.out.println("percent done"+ (double)(newnode-numseq) /(double)numseq * 100);
            a = tree[i];
            b = tree[i+1];
            len_a = si.sl[a];
            len_b = si.sl[b];
            if (len_a < len_b){
                j = a;
                a = b;
                b = j;
                j = len_a;
                len_a = len_b;
                len_b = j;
            }
            dp = reInitializeDPStructure(dp,len_a,len_b);
               //add_poits_to_matrix
               // dp = dp_matrix_init(dp,len_a,len_b);
                add_ptm(matches,dp.m,a,b);
             //   MatrixHelper.printMatrix(dp.m);
                dp = consistency_check(dp,len_a,len_b,dia);
//                //add_ptm2(matches,dp,a,b);
//                //dp = consistency_check2(dp,len_a,len_b,dia);


            path = new int[len_a+len_b+2];
            for (j = len_a+len_b+1;j>=0;j--){
                path[j] = 0;
            }

           // System.out.println(a);
//            int[] temp = profile[a];
//            int[] temp2= si.s[a];

            if (a < numseq){
				System.out.println("sequence a");
				MatrixHelper.printArray(si.s[a]);
				System.out.println(si.s[a].length + " "+ len_a + " " + a);

				profile[a] = make_profile(profile[a],si.s[a],len_a,submatrix);

            }
            if (b < numseq){
                profile[b] = make_profile(profile[b],si.s[b],len_b,submatrix);
            }
			System.out.println("lena = " + len_a + ", lenb = "+len_b);
			set_gap_penalties(profile[a],len_a,si.nsip[b]);
            set_gap_penalties(profile[b],len_b,si.nsip[a]);
			profa = profile[a];
            profb = profile[b];
            path = main_fast_dyn(path,dp,profa,profb,len_a,len_b);
			System.out.println("i = " + i + ", a = " + a + ", b = " + b);
			MatrixHelper.printArray(path);

            profile[newnode] = new int [64*(path[0]+1)];
			System.out.println("newnode " + newnode);
			si = update(si,profile,a,b,newnode,path);
			System.out.println("sl after update");
			MatrixHelper.printArray(si.sl);
//            if (points){
//                update_hash(si,matches,a,b,newnode);
//            }
//            free(profa);
//            free(profb);
            newnode++;
        }
//		if(!quiet)fprintf(stderr,"\r%8.0f percent done",100.0);
//		if(!quiet)fprintf(stderr,"\n");
//		print_alignment(si,outfile);
//		if (!df){
//			for (i = numprofiles;i--;){
//				free(dm[i]);
//			}
//			free(dm);
    }



	static SequenceInfo update(SequenceInfo si,int[][] profile,int a,int b,int newnode,int[] path)
	{
		int i,c;
		int posa = 0;
		int posb = 0;
		int adda = 0;
		int addb = 0;
		int[] gap_a ;
		int[] gap_b ;
		int len_a;
		int len_b;
		int pos_profa=0;
		int pos_profb=0;
		int pos_newp=0;
//	int* newp = 0;
//	int* profa = 0;
//	int* profb = 0;
		len_a = si.sl[a];
		len_b = si.sl[b];

//		newp = profile[newnode];
//		profa = profile[a];
//		profb = profile[b];
		si.sl[newnode] = path[0];
		si.relpos[newnode] = new int[si.sl[newnode]+1];
		for (i = si.sl[newnode];i>=0;i--){
			si.relpos[newnode][i] = i;
		}
		gap_a = new int[len_a+1];
		gap_b = new int[len_b+1];

		for (i = len_a;i>=0;i--){
			gap_a[i] = 0;
		}
		for (i = len_b;i>=0;i--){
			gap_b[i] = 0;
		}

		c = 1;
		while(path[c] != 3){
			if (path[c]==0){
				//		fprintf(stderr,"Align	%d\n",c);
				for (i = 64-1;i>=0; i--){
					profile[newnode][i+pos_newp] = profile[a][i+pos_profa] + profile[b][i+pos_profb];
				}
				si.relpos[a][posa] += adda;
				si.relpos[b][posb] += addb;
				pos_profa += 64;
				pos_profb += 64;
				posa++;
				posb++;
			}
			if ((path[c] & 1)!=0){
				//		fprintf(stderr,"Gap_A:%d\n",c);
				for (i = 64-1;i>=0; i--){
					profile[newnode][i+pos_newp] = profile[b][i+pos_profb];
				}
				si.relpos[b][posb] += addb;
				adda += 1;
				gap_a[posa] += 1;
				pos_profb += 64;
				posb++;
				if ((path[c] & 4)!=0){
					//			fprintf(stderr,"Gap_open");
					profile[newnode][9+pos_newp] += si.nsip[a];//1;
					i = si.nsip[a] *gpo;
					profile[newnode][32+pos_newp] -=  i;//a
					profile[newnode][33+pos_newp] -=  i;//b
					profile[newnode][34+pos_newp] -=  i;//c
					profile[newnode][35+pos_newp] -=  i;//d
					profile[newnode][36+pos_newp] -=  i;//e
					profile[newnode][37+pos_newp] -=  i;//f
					profile[newnode][38+pos_newp] -=  i;//g
					profile[newnode][39+pos_newp] -=  i;//h
					profile[newnode][40+pos_newp] -=  i;//i
					//newp[41] -=  i;//j
					profile[newnode][42+pos_newp] -=  i;//k
					profile[newnode][43+pos_newp] -=  i;//l
					profile[newnode][44+pos_newp] -=  i;//m
					profile[newnode][45+pos_newp] -=  i;//n
					//newp[46] -=  i;//o
					profile[newnode][47+pos_newp] -=  i;//p
					profile[newnode][48+pos_newp] -=  i;//q
					profile[newnode][49+pos_newp] -=  i;//r
					profile[newnode][50+pos_newp] -=  i;//s
					profile[newnode][51+pos_newp] -=  i;//t
					//newp[52] -=  i;//u
					profile[newnode][53+pos_newp] -=  i;//v
					profile[newnode][54+pos_newp] -=  i;//w
					profile[newnode][52+pos_newp] -=  i;//x
					profile[newnode][53+pos_newp] -=  i;//y
					profile[newnode][54+pos_newp] -=  i;//z
				}
				if ((path[c] & 16)!=0){
					profile[newnode][9+pos_newp] += si.nsip[a];//1;
					i = si.nsip[a] *gpo;
					profile[newnode][32+pos_newp] -=  i;//a
					profile[newnode][33+pos_newp] -=  i;//b
					profile[newnode][34+pos_newp] -=  i;//c
					profile[newnode][35+pos_newp] -=  i;//d
					profile[newnode][36+pos_newp] -=  i;//e
					profile[newnode][37+pos_newp] -=  i;//f
					profile[newnode][38+pos_newp] -=  i;//g
					profile[newnode][39+pos_newp] -=  i;//h
					profile[newnode][40+pos_newp] -=  i;//i
					//newp[41] -=  i;//j
					profile[newnode][42+pos_newp] -=  i;//k
					profile[newnode][43+pos_newp] -=  i;//l
					profile[newnode][44+pos_newp] -=  i;//m
					profile[newnode][45+pos_newp] -=  i;//n
					//newp[46] -=  i;//o
					profile[newnode][47+pos_newp] -=  i;//p
					profile[newnode][48+pos_newp] -=  i;//q
					profile[newnode][49+pos_newp] -=  i;//r
					profile[newnode][50+pos_newp] -=  i;//s
					profile[newnode][51+pos_newp] -=  i;//t
					//newp[52] -=  i;//u
					profile[newnode][53+pos_newp] -=  i;//v
					profile[newnode][54+pos_newp] -=  i;//w
					profile[newnode][52+pos_newp] -=  i;//x
					profile[newnode][53+pos_newp] -=  i;//y
					profile[newnode][54+pos_newp] -=  i;//z
				}
			}
			if ((path[c] & 2)!=0){
				//		fprintf(stderr,"Gap_B:%d\n",c);
				System.out.println("size of profile "+ profile.length + " " + profile[newnode].length);
				for (i = 64-1;i>=0; i--){
					profile[newnode][i+pos_newp] = profile[a][i+pos_profa];
				}
				si.relpos[a][posa] += adda;
				addb += 1;
				gap_b[posb] += 1;
				posa++;
				pos_profa+=64;
				if ((path[c] & 4)!=0){
					//			fprintf(stderr,"Gap_open");
					profile[newnode][9+pos_newp] += si.nsip[b];//1;
					i = si.nsip[b] *gpo;
					profile[newnode][32+pos_newp] -=  i;//a
					profile[newnode][33+pos_newp] -=  i;//b
					profile[newnode][34+pos_newp] -=  i;//c
					profile[newnode][35+pos_newp] -=  i;//d
					profile[newnode][36+pos_newp] -=  i;//e
					profile[newnode][37+pos_newp] -=  i;//f
					profile[newnode][38+pos_newp] -=  i;//g
					profile[newnode][39+pos_newp] -=  i;//h
					profile[newnode][40+pos_newp] -=  i;//i
					//newp[41] -=  i;//j
					profile[newnode][42+pos_newp] -=  i;//k
					profile[newnode][43+pos_newp] -=  i;//l
					profile[newnode][44+pos_newp] -=  i;//m
					profile[newnode][45+pos_newp] -=  i;//n
					//newp[46] -=  i;//o
					profile[newnode][47+pos_newp] -=  i;//p
					profile[newnode][48+pos_newp] -=  i;//q
					profile[newnode][49+pos_newp] -=  i;//r
					profile[newnode][50+pos_newp] -=  i;//s
					profile[newnode][51+pos_newp] -=  i;//t
					//newp[52] -=  i;//u
					profile[newnode][53+pos_newp] -=  i;//v
					profile[newnode][54+pos_newp] -=  i;//w
					profile[newnode][52+pos_newp] -=  i;//x
					profile[newnode][53+pos_newp] -=  i;//y
					profile[newnode][54+pos_newp] -=  i;//z
				}
				if ((path[c] & 16)!=0){
					//			fprintf(stderr,"Gap_close");
					profile[newnode][9+pos_newp] += si.nsip[b];//1;
					i = si.nsip[b] *gpo;
					profile[newnode][32+pos_newp] -=  i;//a
					profile[newnode][33+pos_newp] -=  i;//b
					profile[newnode][34+pos_newp] -=  i;//c
					profile[newnode][35+pos_newp] -=  i;//d
					profile[newnode][36+pos_newp] -=  i;//e
					profile[newnode][37+pos_newp] -=  i;//f
					profile[newnode][38+pos_newp] -=  i;//g
					profile[newnode][39+pos_newp] -=  i;//h
					profile[newnode][40+pos_newp] -=  i;//i
					//newp[41] -=  i;//j
					profile[newnode][42+pos_newp] -=  i;//k
					profile[newnode][43+pos_newp] -=  i;//l
					profile[newnode][44+pos_newp] -=  i;//m
					profile[newnode][45+pos_newp] -=  i;//n
					//newp[46] -=  i;//o
					profile[newnode][47+pos_newp] -=  i;//p
					profile[newnode][48+pos_newp] -=  i;//q
					profile[newnode][49+pos_newp] -=  i;//r
					profile[newnode][50+pos_newp] -=  i;//s
					profile[newnode][51+pos_newp] -=  i;//t
					//newp[52] -=  i;//u
					profile[newnode][53+pos_newp] -=  i;//v
					profile[newnode][54+pos_newp] -=  i;//w
					profile[newnode][52+pos_newp] -=  i;//x
					profile[newnode][53+pos_newp] -=  i;//y
					profile[newnode][54+pos_newp] -=  i;//z
				}
			}
			pos_newp += 64;
			c++;
		}
		for (i = 64-1;i>=0; i--){
			profile[newnode][i+pos_newp] = 0;
		}

		//fprintf(stderr,"%d-%d	%d	%d\n",c,path[0],len_a,len_b);
		si.nsip[newnode] = si.nsip[a] + si.nsip[b];
		si.sip[newnode] = new int[si.nsip[newnode]];
		c =0;
		for (i = si.nsip[a]-1;i>=0;i--){
			si.sip[newnode][c] = si.sip[a][i];
			update_gaps(si.sl[si.sip[a][i]],si.gis[si.sip[a][i]],si.sl[newnode],gap_a);
			c++;
		}
		for (i = si.nsip[b]-1;i>=0;i--){
			si.sip[newnode][c] = si.sip[b][i];
			update_gaps(si.sl[si.sip[b][i]],si.gis[si.sip[b][i]],si.sl[newnode],gap_b);
			c++;
		}
	
		return si;
	}

	static void update_gaps(int old_len,int[] gis,int new_len,int[] newgaps)
	{
		int i,j;
		int add = 0;
		int rel_pos = 0;
		for (i = 0; i <= old_len;i++){
			add = 0;
			for (j = rel_pos;j <= rel_pos + gis[i];j++){
				if (newgaps[j] != 0){
					add += newgaps[j];
				}
			}
			rel_pos += gis[i]+1;
			gis[i] += add;
		}
	}

    static void set_gap_penalties(int[] prof,int len,int nsip)
    {
        int i;
        int pos=0;
        pos+=  (64 *(len));
        i = len;
        while(i-->0){
            pos -= 64;
            prof[26+pos] = prof[41+pos]*nsip;
            prof[27+pos] = prof[46+pos]*nsip;
        }
    }

    static int[] main_fast_dyn(int[] path, DPStructure dp,int[] prof1,int[] prof2,int len_a,int len_b)
    {
        int i,j,c;
    	int i_limit = 0;
    	int j_limit = 0;
    	int startx = 0;
    	int starty = 0;
    	int endx = 0;
    	int endy = 0;
    	//int* tx = 0;
    	//int* ty = 0;
    //	int** trace = 0;
    //	int* gap_a = 0;
   // 	int* gap_b = 0;
    //	int* align = 0;
    //	int* tracep = 0;
    	int pa = 0;
    	int pga = 0;
    	int pgb = 0;
    	int ca = 0;
    	int cga = 0;
    	int cgb = 0;
    	int ipos;
    	int jpos;
    	int[] freq = new int[(len_a+1)*26];

    //	tx = dp->true_x;
    //	ty = dp->true_y;


    	//freq = tmalloc((len_a+1) * 26 * sizeof(unsigned int));
    	int pos_prof1 =0;
    	int pos_prof2=0;
    	int pos_freq =0;

    	pos_prof1 +=  len_a << 6;
    	pos_freq += len_a *26;
    	for (i = len_a-1;i>=0;i--){
    		pos_prof1 -= 64;
    		pos_freq -= 26;
    		c = 1;
    		for (j = 26-1;j>=0; j--){
    			if(prof1[pos_prof1+j]!=0){
    				freq[pos_freq+c] = j;
    				c++;
    			}
    		}
    		freq[pos_freq+0] = c;
    	}

    //	align = dp->a;
    //	gap_a = dp->ga;
    //	gap_b = dp->gb;
    	dp.a[0] = 0;
    	dp.ga[0] = -INFINITY;
    	dp.gb[0] = -INFINITY;
    //	trace = dp->tb;
    	endx = len_a;
    	startx = len_a;
    	endy = len_b;
    	starty = len_b;

    	dp.tb[len_a][len_b] = 32;

    	pos_prof1 +=  len_a << 6;

    	pos_freq += len_a *26;

    	do{
    		while(dp.true_x[startx] != 2){
    			startx--;
    		}
    		while(dp.true_y[starty] != 2){
    			starty--;
    		}
    		i_limit = endx-startx;
    		j_limit = endy-starty;
    		//copy last cell;
    		dp.a[j_limit] = dp.a[0];
    		dp.ga[j_limit] = dp.ga[0];
    		dp.gb[j_limit] = dp.gb[0];
    		//init of first row;
    		//tracep = trace[endx];
    		j = j_limit;
    		if (endx == len_a){
    			while(--j>0){
    				dp.a[j] = -INFINITY;
    				dp.ga[j] = 0;
    				dp.gb[j] = -INFINITY;
    				dp.tb[endx][starty+j] = 8;
    			}
    		}else{
    			pos_prof2 +=  endy << 6;
    			while(--j>0){
    				pos_prof2 -= 64;
    				dp.tb[endx][starty+j] = 0;
    				dp.a[j] = -INFINITY;
    				dp.ga[j] = dp.a[j+1] + prof2[pos_prof2+26];
    				if (dp.ga[j+1] > dp.ga[j]){
    					dp.ga[j] = dp.ga[j+1];
    					dp.tb[endx][starty+j] |= 8;
    				}
    				dp.gb[j] = -INFINITY;
    			}
    			pos_prof2 -= (starty+1) << 6;//+1 cos while(--j) stops at 1;(1-1 = 0 stop!!!)
    		}
    		dp.a[0] = -INFINITY;
    		dp.ga[0] = -INFINITY;
    		dp.gb[0] = -INFINITY;
    		pos_prof2 += starty << 6;
    		i = i_limit;
    		while(--i>0){
    			pos_prof1 -= 64;

    			pos_freq -= 26;

    			ipos = startx+i;
    			//tracep = trace[ipos];

    			pa = dp.a[j_limit];
    			pga = dp.ga[j_limit];
    			pgb = dp.gb[j_limit];

    			dp.a[j_limit] = -INFINITY;
    			dp.ga[j_limit] = -INFINITY;

    			dp.tb[ipos][endy] = 0;

    			if (endy == len_b){
    				dp.gb[j_limit] = 0;
    				dp.tb[ipos][endy] |= 16;
    			}else{
    				dp.gb[j_limit] = pa+prof1[pos_prof1+26];
    				if(pgb > dp.gb[j_limit]){
    					dp.gb[j_limit] = pgb;//pgb+prof2[endy][27];
    					dp.tb[ipos][endy] |= 16;
    				}
    			}
    			j = j_limit;
    			pos_prof2 += j_limit << 6;
    			while(--j>0){
    				pos_prof2 -= 64;
    				jpos = starty+j;
    				ca = dp.a[j];
    				cga = dp.ga[j];
    				cgb = dp.gb[j];
    				dp.a[j] = pa;
    				dp.tb[ipos][jpos] = 1;
    				if((c = pga+prof2[pos_prof2+ 64+26]) > dp.a[j]){ // TODO: check prof2 value
    					dp.a[j] = c;//pga+prof1[ipos+1][26];
    					dp.tb[ipos][jpos] = 2;
    				}
    				if((c = pgb+prof1[pos_prof1+ 64+26]) > dp.a[j]){
    					dp.a[j] = c;//pgb+prof2[jpos+1][26];
    					dp.tb[ipos][jpos] = 4;
    				}
    				for (c = freq[pos_freq+0]-1;c>0;c--){
    					dp.a[j] += prof1[pos_prof1+freq[pos_freq+c]]*prof2[pos_prof2+freq[pos_freq+c] | 32];
    				}
    				dp.ga[j] = dp.a[j+1]+prof2[pos_prof2+26];
    				if (dp.ga[j+1] > dp.ga[j]){
    					dp.ga[j] = dp.ga[j+1];//gap_a[j+1]+prof1[ipos][27];
    					dp.tb[ipos][jpos] |= 8;
    				}
    				dp.gb[j] = ca+prof1[pos_prof1+26];// prof2[jpos][26];
    				if(cgb > dp.gb[j]){
    					dp.gb[j] = cgb;//cgb+prof2[jpos][27];
    					dp.tb[ipos][jpos] |= 16;
    				}
    				pa = ca;
    				pga = cga;
    				pgb = cgb;
    			}
    			pos_prof2 -= 64;
    			//LAST CELL (0)
    			ca = dp.a[0];
    			cga = dp.ga[0];
    			cgb = dp.gb[0];

    			dp.a[0] = pa;
    			dp.tb[ipos][starty] = 1;
    			if((c = pga+prof2[pos_prof2+ 64+26]) > dp.a[0]){
    				dp.a[0] = c;//pga+prof1[ipos+1][26];
    				dp.tb[ipos][starty] = 2;
    			}
    			if((c = pgb+prof1[pos_prof1+ 64+26]) > dp.a[0]){
    				dp.a[0] = c;//pgb+prof2[jpos+1][26];
    				dp.tb[ipos][starty] = 4;
    			}
    			for (c = freq[pos_freq+0]-1;c>0;c--){
    				dp.a[j] += prof1[pos_prof1+freq[pos_freq+c]]*prof2[pos_prof2+freq[pos_freq+c] | 32];
    			}

    			dp.ga[j] = -INFINITY;

    			dp.gb[0] = ca+prof1[pos_prof1+26];
     			if(cgb > dp.gb[0]){
    				dp.gb[0] = cgb;
    				dp.tb[ipos][starty] |= 16;
    			}
    		}
    		pos_prof1 -= 64;

    		pos_freq -= 26;
    		//tracep = trace[startx];
    		j = j_limit;

    		pos_prof2 += j_limit << 6;
    		pa = dp.a[j];
    		pga = dp.ga[j];
    		pgb = dp.gb[j];

    		dp.a[j] = -INFINITY;
    		dp.ga[j] = -INFINITY;
    		dp.gb[j_limit] = -INFINITY;
    		while(--j>0){
    			pos_prof2 -= 64;

    			jpos = starty+j;

    			ca = dp.a[j];
    			cga = dp.ga[j];
    			cgb = dp.gb[j];

    			dp.a[j] = pa;
    			dp.tb[startx][jpos] = 1;
    			if((c = pga+prof2[pos_prof2+ 64+26]) > dp.a[j]){
    				dp.a[j] = c;//pga+prof1[ipos+1][26];
    				dp.tb[startx][jpos] = 2;
    			}
    			//Gap_b->Align
    			if((c = pgb+prof1[pos_prof1+ 64+26]) > dp.a[j]){
    				dp.a[j] = c;//pgb+prof2[jpos+1][26];
    				dp.tb[startx][jpos] = 4;
    			}

    			for (c = freq[pos_freq+0]-1;c>0;c--){
    				dp.a[j] += prof1[pos_prof1 + freq[pos_freq+c]]*prof2[pos_prof2+freq[pos_freq+c] | 32];
    			}
    			dp.ga[j] = dp.a[j+1]+prof2[pos_prof2+26];
    			if (dp.ga[j+1] > dp.ga[j]){
    				dp.ga[j] = dp.ga[j+1];//gap_a[j+1]+prof1[ipos][27];
    				dp.tb[startx][jpos] |= 8;
    			}
    			dp.gb[j] = -INFINITY;
    			pa = ca;
    			pga = cga;
    			pgb = cgb;
    		}

    		pos_prof2 -= 64;

    		ca = dp.a[0];
    		cga = dp.ga[0];
    		cgb = dp.gb[0];
    		dp.a[0] = pa;
    		dp.tb[startx][starty] = 1;
    		if((c = pga+prof2[pos_prof2+ 64+26]) > dp.a[0]){
    			dp.a[0] = c;//pga+prof1[ipos+1][26];
    			dp.tb[startx][starty] = 2;
    		}
    		if((c = pgb+prof1[pos_prof1+ 64+26]) > dp.a[0]){
    			dp.a[0] = c;//pgb+prof2[jpos+1][26];
    			dp.tb[startx][starty] = 4;
    		}

    		for (c = freq[pos_freq+0]-1;c>0;c--){
    			dp.a[j] += prof1[pos_prof1+freq[pos_freq+c]]*prof2[pos_prof2+freq[pos_freq+c] | 32];
    		}
    		dp.ga[j] = dp.a[j+1]+prof2[pos_prof2+26];
    		//fprintf(stderr,"Gap-a:%d\n",prof2[26]);
    		//Gap_a->Gap_a
    		if (dp.ga[j+1] > dp.ga[j]){
    			dp.ga[j] = dp.ga[j+1];//gap_a[j+1]+prof1[ipos][27];
    			dp.tb[startx][starty] |= 8;
    		}
    		dp.gb[0] = ca+prof1[pos_prof1+26];// prof2[jpos][26];
    		//fprintf(stderr,"Gap-b:%d\n",prof1[26]);
    		//Gap_b->Gap_b
    		if(cgb > dp.gb[0]){
    			dp.gb[0]= cgb;
    			dp.tb[startx][starty] |= 16;
    		}
    		pos_prof2 -= starty<<6;
    		//fprintf(stderr,"\nMOVED-::%d\n",(starty) << 6);
    		endx = startx;
    		endy = starty;
    		startx--;
    		starty--;
    	}while (startx >= 0 || starty >= 0);

    	//free(freq);

    	ca = dp.gb[0];
    	c = 2;
    	if(dp.ga[0] > ca){
    		ca = dp.ga[0];
    		c = 1;
    	}
    	if(dp.a[0] >= ca){
    		//ca = align[0];
    		c = 0;
    	}
    	//fprintf(stderr,"STATE:%d	%d\n",c,ca);
    	ca = c;

    	i = 0;
    	j = 0;
    	c = 1;
    	while(dp.tb[i][j] < 32){
    	//	fprintf(stderr,"%d->%d	%d:%d	%d:%d\n",c,trace[i][j],i,j,len_a,len_b);
    		switch(ca){
    			case 0:
    				if ((dp.tb[i][j] & 2)!=0){
    					ca = 1;
    					if(i+1!= len_a){
    						path[c+1] |= 16;
    	//					fprintf(stderr,"GAP_CLOSE\n");
    					}
    				}else if ((dp.tb[i][j] & 4)!=0){
    					ca = 2;
    					if(j+1!= len_b){
    						path[c+1] |= 16;
    	//					fprintf(stderr,"GAP_CLOSE\n");
    					}
    				}

    				//path[c] = 0;
    				i++;
    				j++;
    			break;
    			case 1:
    				if((dp.tb[i][j] & 8)!=0){
    					ca = 1;
    					if(i!=0 && i!= len_a){
    	//					fprintf(stderr,"GAP_EXT\n");
    						if((path[c]&16)==0){ //TODO: check value
    							path[c] |= 8;
    						}
    					}
    				}else{
    					ca = 0;
    					if(i!=0 && i!= len_a){
    	//					fprintf(stderr,"GAP_OPEN\n");
    						path[c] |= 4;
    					}

    				}
    				path[c] |= 1;
    				j++;
    			break;
    			case  2:
    				if((dp.tb[i][j] & 16)!=0){
    					ca = 2;
    					if(j !=0 && j != len_b){
    	//					fprintf(stderr,"GAP_EXT\n");
    						if((path[c]&16)==0){
    							path[c] |= 8;
    						}
    					}
    				}else{
    					ca = 0;
    					if(j!=0 && j != len_b){
    	//					fprintf(stderr,"GAP_OPEN\n");
    						path[c] |= 4;
    					}
    				}
    				path[c] |= 2;
    				i++;
    			break;
    		}
    		c++;
    		//System.out.println("i is "+i);
    		//System.out.println("j is "+j);
    	}
    	path[0] = c-1;
    	path[c] = 3;
    	return path;
    }

    static int[] make_profile(int[] prof,int[] seq,int len,int[][] subm)
    {
        int i,j,c;
        prof = new int[(len+1)*64];
        int pos =0;
        pos+=  (64 *len);
        //fprintf(stderr,"Len:%d	%d\n",len,64*len);
        for (i = 64-1;i>=0;i--){
            prof[pos+i] = 0;
        }
        prof[pos + 9+32] = subm[0][9];
        prof[pos + 14+32] = subm[0][14];
        i = len;
       // pos = 64*len;
        while(i-->0){
            pos -= 64;
            //fprintf(stderr,"-64\n");
            for (j = 26-1;j>=0; j--){
                prof[pos+j] = 0;
            }
            c = seq[i];
            prof[pos+c] +=1;
            pos += 32;
            for(j = 32-1;j>=0;j--){
                prof[pos+j] = subm[c][j];
            }
            pos -= 32;
        }
        return prof;
    }

    static DPStructure consistency_check(DPStructure dp,int len_a,int len_b, int dia)
    {

        int i,j;
        int c = 0;


       // tbop = tb[len_a];
        for (j = len_b;j>=0;j--){
            dp.tb[len_a][j] = 0;
        }
        //tbop[len_b] = 0;
        for (i = len_a-1;i>=0;i--){
//            mxp = mx[i];
//            tbp = tb[i];
            dp.tb[i][len_b] = 0;
            for (j = len_b-1;j>=0;j--){
                dp.tb[i][j] = 0;
                if (dp.m[i][j]!=0){
                    dp.tb[i][j] = dp.tb[i+1][j+1] + 1;
                }
                dp.m[i][j] += dp.m[i+1][j+1];
                if (dp.m[i][j+1] > dp.m[i][j]){
                    dp.m[i][j] = dp.m[i][j+1];
                    dp.tb[i][j] = -1;
                }
                if (dp.m[i+1][j] > dp.m[i][j]){
                    dp.m[i][j] = dp.m[i+1][j];
                    dp.tb[i][j] = -2;
                }
            }
//            mxop = mxp;
//            tbop = tbp;
        }
        c = 0;
        i = 0;
        j = 0;
        while (i < len_a && j < len_b){
            //fprintf(stderr,"%d	%d\n",tb[i][j] & 0x0000ffff,c);
            switch (dp.tb[i][j]){
                case -1:
                    //printf("%d:%d	%d\n",i,j,c);
                    c = 0;
                    j++;
                    break;
                case -2:
                    //printf("%d:%d	%d\n",i,j,c);
                    c = 0;
                    i++;
                    break;
                default:
                    if (c!=0){
                        c = dp.tb[i][j];
                        if (c < dia){
                            c = 0;
                        }else{
                            c -= 1;
                            while (--c>0){
                                dp.true_x[i+c] = 2;
                                dp.true_y[j+c] = 2;
                            }
                        }
                    }
                    //	tx[i] = 2;
                    //	ty[j] = 2;
                    i++;
                    j++;
                    break;

            }
        }
        //exit(0);
        dp.true_x[0] = 2;
        dp.true_y[0] = 2;

        return dp;
    }

    static void add_ptm(int[][][] matches,int[][] matrix,int a,int b)
    {
        int i,j,c;
        int[] posa = null;
        int[] posb = null;
         for (c =8000-1; c>=0;c--){
            if (matches[c]!=null){
                if (matches[c][a] != null && matches[c][a][0]!= 1){
                    if (matches[c][b] !=null &&matches[c][b][0]!= 1){
                        posa = matches[c][a];
                        posb = matches[c][b];
                       // System.out.println(matrix[]);
                        for (i = posa[0]-1;i>0;i--){
                            for (j = posb[0]-1;j>0;j--){
                                //System.out.println("pos a " + posa[i] + "posb " + posb[j] );
                             //   System.out.println("size of matrix " + matrix[posa[i]] .length);
                                matrix[posa[i]][posb[j]] += 1;
                            }
                        }
                    }
                }
            }
        }
    }

    static DPStructure dp_matrix_init(DPStructure dp,int x,int y)
    {
        int[][] mx = null;
            //int** tb = 0;
        int[] mxp = null;
            //int* tbp = 0;
        int[] tx = null;
        int[] ty = null;
        int i,j;

        tx = dp.true_x;
        ty = dp.true_y;
        //tb = dp->tb;
        mx = dp.m;
        for (i = x;i>=0;i--){
            //tbp = tb[i];
            mxp = mx[i];
            tx[i] = 0;
            for (j = y;j>=0;j--){
                //	tbp[j] = 0;
                mxp[j] = 0;
            }
        }
        for (j = y;j>=0;j--){
            ty[j] = 0;
        }
        return dp;
    }
   static DPStructure reInitializeDPStructure(DPStructure dp,int x,int y)
    {
        int i;
        if ( x > dp.x || y > dp.y){
            i = 1;
            while (i <= y){
                i <<= 1;
                //	printf("i:%d	y:%d\n",i,y);
            }
            y = i-1;
            i = 1;
            while (i <= x){
                i <<= 1;
                //printf("i:%d	y:%d\n",i,y);
            }
            x = i-1;
            //printf("NEWX:%d	NEWY:%d\n",x,y);
            dp.a = new int[y+1];
            dp.ga = new int[y+1];
            dp.gb = new int[y+1];
            dp.tb = new int[x+1][y+1];
            dp.m = new int [x+1][y+1];
//            for ( i = 1; i <= x;i++) {
//                dp.tb[i] =  new int [1];
//                int[] temp1 = new int[1];
//                temp1[0] = (i * (y + 1));
//                dp.tb[i] = temp1;
//                dp.m[i] = temp1;
//            }

            dp.true_x = new int[x+1];
            dp.true_y = new int[y+1];
            dp.x = x;
            dp.y = y;
        }
        return dp;
    }

    static DPStructure initializeDPStructure(DPStructure dp, int x, int y) {

        int i;
        dp.x = x;
        dp.y = y;
        dp.a = new int[y+1];
        dp.ga = new int[y+1];
        dp.gb = new int[y+1];
        dp.true_x = new int[x+1];
        dp.true_y = new int[y+1];
        dp.tb = new int[x+1][y+1];
        dp.m = new int [x+1][y+1];

       // TODO : Check how to move the pointer index
//        for ( i = 1; i <= x;i++) {
//           //  dp.tb[i] =  new int [1];
//            int[] temp1 = new int[1];
//            temp1[0] = (i * (y + 1));
//            dp.tb[i] = temp1;
//            dp.m[i] = temp1;
//            System.out.println(dp.tb[i][0]);
//        }

        return dp;
    }

    static int[][] populateScoringMatrix(int[] scoringArray,int[][] scoringmatrix,int gpo,char[] letters)
    {
        int i,j;
        int m_pos = 0;
        m_pos = 0;
        for (i = 0;i < 23;i++){
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
        for (i = 26-1;i>=0;i--){
            scoringmatrix[i][9] = -gpo;
        }
        return scoringmatrix;
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
            max = -INFINITY;
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
                       // System.out.println("Printing p's");
                        for (i = (p1[0] & 0x0000ffff)-1;i>0;i--){
                           // int temp_p = p[temp+p1[i]];
                           // p = diagonal - p1[i];
                            for (j = (p2[0] & 0x0000ffff)-1;j>0;j--){
                                diagonal[len_a-p1[i] + p2[j]] += tmp2;
                               // System.out.print( diagonal[len_a-p1[i] + p2[j]] + " ");
                            }
                           // System.out.println("\n");
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
        //    MatrixHelper.printArray(si.s[i]);
            si.sl[i] = key.length();
            si.sn[i] = value;
            si.lsn[i] = value.length();
            i++;
        }

        si.sl[numseq] = si.sl[numseq-1];
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
            si.sip[i] = new int[1];
            si.nsip[i] = 1;
            si.sip[i][0] = i;
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
     //   System.out.println("Printing Keys\n");
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
                //    System.out.println("key " + key + " i " + i + " j " + j);
                    matches[key][i] = new int[10];
                    matches[key][i][0] = 1;
                }
                if ((matches[key][i][0] %10) == 0){
                 //  System.out.println("Realloc key " + key + " i " + i + " j " + j);
                    matches[key][i] = Arrays.copyOf(matches[key][i], (matches[key][i][0] + 10));
                }
                matches[key][i][matches[key][i][0]] = j;
                matches[key][i][0] += 1;
            }
            //System.out.println("\n");
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
                      //  System.out.println("c = 10 ");
                      //  MatrixHelper.printMatrix(out);
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
//                System.out.println("c!=0");
//                MatrixHelper.printMatrix(out);
                for (f = 0;f < c;f++){
                    matches[map[f]][i] = out[f];
                    matches[map[f]][i][0] |= 0x00010000;
                }
                //free(out);
                c = 0;
            }
        }
       // MatrixHelper.printArray(patterns);

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


