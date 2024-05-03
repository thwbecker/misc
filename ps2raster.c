/*--------------------------------------------------------------------
 *
 * ps2raster Converts one or several ps file(s) to raster formats using GhostScript
 * It works by modifying the BoundingBox and PageSize in order that the raster
 * image will have a size are determined by the BoundingBox. 
 * As an option, a tight BoundingBox may be computed.
 * ps2raster uses the ideas of the EPS2XXX.m from Primoz Cermelj published in MatLab Central
 * and of psbox_sh of Remko Scharroo
 *
 * To compile, run:
 *	gcc ps2raster.c -lm -o ps2raster	(in *nix like environments)
 *      cl ps2raster.c				(with VC)
 *
 *--------------------------------------------------------------------*/
/*
 * Author:	Joaquim Luis
 * Created:	15-FEB-2005
 * Modified:	15-NOV-2005 (It used to fail too many times in detecting if is GMT_ps)
 *		
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define	FALSE	0
#define	TRUE	1


void run_cmd(char *name, char *cmd, int rem1, int rem2);

main(int argc, char **argv) {
	int i, j, k, len, x0, x1, y0, y1, w, h, n_files = 0, n_alloc = 100;
	int error = FALSE, adjust_BB = FALSE, no_rem = FALSE, have_list = FALSE, verbose = FALSE;
	int ps_GMT, got_BB, got_adjusted_BB, atend;
	
	char ps_file[BUFSIZ];
	char **ps_names;
	char line[BUFSIZ], c[BUFSIZ], c1[10], c2[10], c3[10], c4[10], c5[10], cmd[BUFSIZ], cmdBB[BUFSIZ];
	char *tmp_file = NULL, *out_file = NULL, *BB_file = NULL, *fullGsPath = NULL;
	char *script_name = NULL, *scriptBB = NULL, *gs_params = NULL, *gs_BB = NULL;
	char *dpi = "300";
	char *ext = ".png";
	char *device = "png16m";
	FILE *fp = NULL, *fpo = NULL, *fpb = NULL, *fpl = NULL;
	
	/* Parameters for all the formats available */
	gs_params = "-q -dNOPAUSE -dBATCH -dDOINTERPOLATE -dUseFlateCompression=true -sDEVICE=";
	gs_BB = "-q -dNOPAUSE -dBATCH -sDEVICE=bbox ";

#ifdef _WIN32
	fullGsPath = "gswin32c ";
	script_name = "ps2raster.bat";
	scriptBB = "to_BB.bat";
#else
	fullGsPath = "gs ";
	script_name = "ps2raster.sc";
	scriptBB = "to_BB.sc";
#endif
	
	/* Check and interpret the command line arguments */
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {
				case '\0':
					break;
				
				/* Supplemental parameters */
				case 'A':
					adjust_BB = TRUE;
					break;
				case 'G':
					fullGsPath = &argv[i][2];
					strcat(fullGsPath," ");		/* We need this blank space */
					break;
				case 'E':
                                        dpi = &argv[i][2];
					break;
				case 'L':
					if ((fpl = fopen (&argv[i][2], "r")) == NULL) {
						fprintf (stderr, "PS2RASTER: Cannot to open list file %s\n", argv[i]); 
						error = TRUE;
					}
					have_list = TRUE;
					break;
				case 'N':
					no_rem = TRUE;
					break;
				case 'T':	/* Select output format */
					for (j = 2; argv[i][j]; j++) {
						switch (argv[i][j]) {
							case 'f':	/* PDF */
								ext = ".pdf";
								device = "pdfwrite";
								break;
							case 'j':	/* JPEG */
								ext = ".jpg";
								device = "jpeg";
								break;
							case 'g':	/* PNG */
								ext = ".png";
								device = "png16m";
								break;
							case 'm':	/* PPM */
								ext = ".ppm";
								device = "ppm16m";
								break;
							case 't':	/* TIFF */
								ext = ".tif";
								device = "tifflzw";
								break;
							default:
								fprintf (stderr, "PS2RASTER ERROR: unrecognized option %s\n", argv[i]); 
								error = TRUE;
								break;
						}
					}
					break;
				case 'V':
					verbose = TRUE;
					break;
				/* Options not recognized */
				default:
					error = TRUE;
					break;
			}
		}
		else {
			strcpy (ps_file, argv[i]);
			n_files++;
		}
	}
	
	if (argc == 1 || error) {
		fprintf (stderr,"ps2raster - Converts one or several PS file(s) to raster formats using GhostScript.\n");
		fprintf (stderr,"It works by modifying the BoundingBox and PageSize in order that the raster\n");
		fprintf (stderr,"image will have a size are determined by the BoundingBox.\n");
		fprintf (stderr,"As an option, a tight BoundingBox may be computed.\n\n");
		fprintf(stderr,"usage: ps2raster <psfile1> <psfile2> <...> [-A] [-E<resolution>]\n");
		fprintf(stderr,"       [-G<ghost_path>] [-L<listfile>] [-N] [-Tf|j|g|m|t] [-V]\n\n");
		
		fprintf (stderr,"	<psfile(s)> postscript file(s) to be converted.\n");
		fprintf(stderr,"\n\tOPTIONS:\n");
		fprintf (stderr,"	-A adjust the BoundingBox to the minimum required by the image contents\n");
		fprintf (stderr,"	-E set raster resolution in dpi [default = 300]\n");
		fprintf (stderr,"	-G Full path to your ghostscript executable.\n");
		fprintf (stderr,"	   NOTE: while in the *nix systems this is generally not necessary,\n");
		fprintf (stderr,"	   under windows ghostscript is not added to the system's path. \n");
		fprintf (stderr,"	   So either you do it yourself, or give the full path here.\n");
		fprintf (stderr,"	   (e.g. -Gc:\\programas\\gs\\gs7.05\\bin\\gswin32c).\n");
		fprintf (stderr,"	-L listfile is an ASCII file with names of ps files to be converted\n");
		fprintf (stderr,"	-N do not remove auxiliary files [default, does]. Auxiliary files are\n");
		fprintf (stderr,"	   build using as base the input ps files and consist of:\n");
		fprintf (stderr,"	   psfile_tmpxx.eps -> postcript with a modified BoundingBox\n");
		fprintf (stderr,"	   ps2raster.sc|bat -> script (or batch) with the GS command that does the job\n");
		fprintf (stderr,"	   Use this option to save the script and run it later with different settings\n");
		fprintf (stderr,"	-T set output raster format [default is png]\n");
		fprintf (stderr,"	   f means PDF\n");
		fprintf (stderr,"	   j means JPEG\n");
		fprintf (stderr,"	   g means PNG\n");
		fprintf (stderr,"	   m means PPM\n");
		fprintf (stderr,"	   t means TIF\n");
		exit (0);
	}

	if (n_files > 1 && have_list) {
		fprintf (stderr, "ps2raster error: Cannot handle both a file list and multiple ps files in input\n");
		error = TRUE;
	}

	if (error) exit (0);

	if (have_list) {		/* Multiple files in a file with their names */
		ps_names = (char **) calloc ((size_t)(n_alloc), sizeof (char *));
		while (fgets (line, BUFSIZ, fpl) != NULL) {
			ps_names[n_files] = (char *) calloc (BUFSIZ, sizeof (char));
			if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
			sscanf (line, "%s",ps_names[n_files]);
			n_files++;
			if (n_files > n_alloc) {
				n_alloc += 100;
				ps_names = realloc (ps_names, n_alloc * sizeof (char *));
			}
		}
		fclose (fpl);
	}
	else if (n_files > 1) {		/* Multiple files given in command line */
		ps_names = (char **) calloc ((size_t)(n_alloc), sizeof (char *));
		j = 0;
		for (k = 1; k < argc; k++) {
			if (argv[k][0] == '-') continue;
			ps_names[j] = (char *) calloc (BUFSIZ, sizeof (char));
			ps_names[j] = argv[k];
			j++;
			if (n_files > n_alloc) {
				n_alloc += 100;
				ps_names = realloc (ps_names, n_alloc * sizeof (char *));
			}
		}
	}
	else {				/* Single file */
		ps_names = (char **) calloc ((size_t)(1), sizeof (char *));
		ps_names[0] = ps_file;
	}

	for (k = 0; k < n_files; k++) {		/* Loop over all input files */
		strcpy(ps_file,ps_names[k]);
		if ((fp = fopen (ps_file, "r")) == NULL) {
			fprintf (stderr, "Cannot to open file %s\n", ps_file);
			continue;
		}	
		ps_GMT = FALSE;
		atend = FALSE;
		got_BB = FALSE;
		got_adjusted_BB = FALSE;

		len = strlen(ps_file);
		j = len - 1;
		for (i = 0; i < len; i++, j--) {
			if (ps_file[j] == '.') break;
		} 
		tmp_file = (char *) calloc ((size_t)(j+12), sizeof(char));
		out_file = (char *) calloc ((size_t)(j+6), sizeof(char));

		/* Adjust to a tight BoundingBox if User required so */
		if (adjust_BB) {
			BB_file = (char *) calloc ((size_t)(j+7), sizeof(char));
			for (i = 0; i < j; i++) BB_file[i] = ps_file[i];
			strcat(BB_file,".bbox");
			sprintf(cmd,"%s",fullGsPath);
			strcat(cmd,gs_BB);
			strcat(cmd,ps_file);
			strcat(cmd," 2> ");
			strcat(cmd,BB_file);
			run_cmd(scriptBB, cmd, 0, 0);	/* Execute the script that computes the tight BB */
			if ((fpb = fopen (BB_file, "r")) == NULL) {
				fprintf (stderr, "Unable to open file %s\n", BB_file);
				exit (0);
			}
			while (fgets (line, BUFSIZ, fpb) != NULL) {
				sscanf (line, "%s",c); 
				if (strcmp(c,"%%BoundingBox:") == 0) {
					sscanf (line, "%s %s %s %s %s",c1,c2,c3,c4,c5); 
					x0 = atoi (c2);		y0 = atoi (c3);
					x1 = atoi (c4);		y1 = atoi (c5);
					got_adjusted_BB = TRUE;
				}
			}
			fclose (fpb);
			run_cmd(BB_file, cmd, 1, 0);	/* Remove the file with BB info */
			free ((void *)BB_file);
		}

		for (i = 0; i < j; i++) tmp_file[i] = ps_file[i];
		for (i = 0; i < j; i++) out_file[i] = ps_file[i];
		strcat(tmp_file,"_tmpxx.eps");
		if ((fpo = fopen (tmp_file, "w")) == NULL) {
			fprintf (stderr, " Unable to open file %s for writing\n", tmp_file);
			exit (0);
		}
		strcat(out_file,ext);

		/* Build the converting ghost command */
		sprintf(cmd,"%s",fullGsPath);
		strcat(cmd,gs_params);
		strcat(cmd,device);
		strcat(cmd," -r");
		strcat(cmd,dpi);
		strcat(cmd," -sOutputFile=");
		strcat(cmd,out_file);
		strcat(cmd," -f ");
		strcat(cmd,tmp_file);

		/* The following assumes that the BoundingBox command will be found in the first 20 lines.
		   It serves for three purposes:
		   1 - See if we have a case of BoundingBox == (atend)
		   2 - Find out if we dealing (or not) with a GMT ps file
		   3 - Speed up the greping code */
		for (i = 0; i < 20; i++) {
			fgets (line, BUFSIZ, fp);
			sscanf (line, "%s",c); 
			if (!got_BB && !strncmp(c,"%%BoundingBox:", 14)) {
				sscanf (line, "%s %s %s %s %s",c1,c2,c3,c4,c5); 
				if (!strncmp(c2,"(atend)", 7)) {
					if (!adjust_BB) {
						fprintf (stderr, "To convert the %s file you must use the -A option\n", ps_file);
						exit(-1); 
					}
					atend = TRUE;
				}
				if (!got_adjusted_BB) {		/* Otherwise BoundingBox was already adjusted above as required by -A option */
					x0 = atoi (c2);		y0 = atoi (c3);
					x1 = atoi (c4);		y1 = atoi (c5);
				}
				w = fabs(x1 - x0);
				h = fabs(y1 - y0);
				fprintf (fpo, "%s 0 0 %d %d\n", c1,w,h);
				got_BB = TRUE;
			}
			else if (!strcmp(line,"%%Creator: GMT\n")) {	/* Not sure that this test allways works in finding GMT ps */
				ps_GMT = TRUE;
				fprintf (fpo, "%s", line);
			}
			else if (!strncmp(c,"%%BeginProlog", 13)) {
				if (!ps_GMT) {		/* It is not a GMT ps file. Lets hope the following will work */
					fprintf (fpo, "<< /PageSize [ %d %d ] >> setpagedevice\n", w,h);
					fprintf (fpo, "gsave %d %d translate\n\n", -x0,-y0);
					fprintf (fpo, "%s", line);
				}
				else
					fprintf (fpo, "%s", line);
			}
			else
				fprintf (fpo, "%s", line);
		}
		if (!got_BB) {
			fprintf (stderr, "Fatal ERROR: the file %s has no BoundingBox in the firts 20 lines\n", ps_file);
			exit(-1); 
		}

		/* Do the rest of the postscript file */
		while (fgets (line, BUFSIZ, fp) != NULL) {
			if (line[0] != '%') {
				fprintf (fpo, "%s", line);
				continue;	/* We are looking only for information on comment lines */
			}
			if (atend) {		/* We still have one more BoundingBox somewhere */
				sscanf (line, "%s",c); 
				if (!strncmp(c,"%%BoundingBox:", 14)) continue;
			}
			if (ps_GMT) {
				if (!strcmp(line,"% Init coordinate system and scales\n")) {
					fprintf (fpo, "<< /PageSize [ %d %d ] >> setpagedevice\n", w,h);
					fprintf (fpo, "gsave %d %d translate\n\n", -x0,-y0);
					fprintf (fpo, "%s", line);
				}
				else
					fprintf (fpo, "%s", line);
			}
			else {
				fprintf (fpo, "%s", line);
			}

		}

		fclose (fpo);
		fclose (fp);

		run_cmd(script_name, cmd, 0, no_rem);	/* Execute the script that converts the PS to raster */
		if (!no_rem)
			run_cmd(tmp_file, cmd, 1, 0);

		free ((void *)tmp_file);
		free ((void *)out_file);
	}

	free ((void *)ps_names);
}

void run_cmd(char *name, char *cmd, int rem1, int rem2) {
	FILE *fp;

	if (rem1) {
#ifdef _WIN32
		sprintf (cmd, "del %s\n", name);
#else
		sprintf (cmd, "rm -f %s\n", name);
#endif
		system(cmd);
		return;
	}

	fp = fopen (name, "w");
#ifdef _WIN32
	fprintf (fp, "echo off\n");
#endif
	fprintf (fp, "%s\n", cmd);
	fclose (fp);

#ifdef _WIN32
	sprintf (cmd, "%s",name);
	system(cmd);
#else
	sprintf (cmd, "sh %s", name);
	system(cmd);
#endif
	if (!rem2) {
#ifdef _WIN32
		sprintf (cmd, "del %s\n", name);
#else
		sprintf (cmd, "rm -f %s\n", name);
#endif
		system(cmd);
	}
}
