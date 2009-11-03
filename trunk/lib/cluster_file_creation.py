import string
import os, sys
import numpy

def CreateXYZ(data, cell, molecules, dir, filename_base, verb=2):
	file = '%s%s%s.xyz' % (dir, os.sep, filename_base)
	try:
		foutput = open(file, 'w')
	except:
		if verb>0:
			print "[ERROR] Could not create %s. Aborting..." % (file)
		sys.exit(1)
		
	tmp = '%d %d\n' % (data.n_frame, len(molecules))
	for ii in molecules:
		tmp += '%d ' % (data.n_atom[ii])
	tmp += '\n'	
	foutput.writelines(tmp)
	tmp = ''
	
	for i in xrange(data.n_frame):
		tmp += 'frame %d\n' % (i)
		for ii in molecules:
			tmp += 'molecule %d\n' % (ii)
			for iii in xrange(data.n_atom[ii]):
				tmp += '%4s %5d %10f %5d %12f %12f %12f\n'\
				% (data.symbol[i][ii][iii], data.atomic_number[i, ii, iii], data.atomic_mass[i, ii, iii], data.atomic_valence[i, ii, iii], data.x[i, ii, iii], data.y[i, ii, iii], data.z[i, ii, iii]) 
		if ( i%100 == 0 ):
			foutput.writelines(tmp)
			tmp = ''
			if verb > 3:
				print "[INFO] Writing frame %d in the file %s" % (i, file)
	
	foutput.write(tmp)
	foutput.close()

def CreateCELL(data, cell, pbc_number, cutoff, dir, filename_base):
	tmp = 'PBC %s' % (pbc_number)
	tmp += 'cutoff %f\n' % (cutoff)
	for i in xrange(data.n_frame):
		tmp += 'frame %d\n' % (i)
		tmp += '%f %f %f %f %f %f\n' % (cell.a[i], cell.b[i], cell.c[i], cell.alpha_deg[i], cell.beta_deg[i], cell.gamma_deg[i])
		tmp += '%.15f %.15f %.15f %.15f %.15f %.15f %.15f\n' % (cell.temp_alpha_cos[i], cell.temp_beta_sin[i], cell.temp_beta_cos[i], cell.temp_gamma_sin[i], cell.temp_gamma_cos[i], cell.temp_beta_term[i], cell.temp_gamma_term[i])

	file = '%s%s%s.cell' % (dir, os.sep, filename_base)
	try:
		foutput = open(file, 'w')
	except:
		print "[ERROR] Could not create %s. Aborting..." % (file)
		sys.exit(1)
	
	foutput.write(tmp)
	foutput.close()

def CreateCM(data, molecules, molecules_J, dir, filename_base):
	tmp = ''
	for i in xrange(data.n_frame):
		tmp += 'frame %d\n' % (i)
		for ii in molecules:
			if ii in molecules_J:
				a = 1
			else:
				a = 0
			tmp += 'molecule %d %d ' % (ii, data.n_electrons[i, ii])
			tmp += '%.15f %.15f %.15f %d\n' % (data.CM_x[i, ii], data.CM_y[i, ii], data.CM_y[i, ii], a)

	file = '%s%s%s.cm' % (dir, os.sep, filename_base)
	try:
		foutput = open(file, 'w')
	except:
		print "[ERROR] Could not create %s. Aborting..." % (file)
		sys.exit(1)
	
	foutput.write(tmp)
	foutput.close()
	
def ScriptFileCreation(dir_cluster, input_dir_cluster, output_dir_cluster, scratch_dir_cluster, zindo_dir_cluster, location_cluster, project_name):
	tmp = ''
	tmp += '#!/bin/bash\n\n'
	tmp += 'DIR="%s"\n' % (dir_cluster)
	tmp += 'INPUT_DIR="%s"\n' % (input_dir_cluster)
	tmp += 'OUTPUT_DIR="%s"\n' % (output_dir_cluster)
	tmp += 'SCRATCH_DIR="%s"\n' % (scratch_dir_cluster)
	tmp += 'ZINDO_DIR="%s"\n\n' % (zindo_dir_cluster)

	tmp += 'if [[ -d $DIR ]]; then\n'
	tmp += '	cd $DIR\n'
	tmp += 'else\n'
	tmp += '	echo "The folder $DIR does not exist, but it is supposed to be the project directory. Aborting..."\n'
	tmp += '	exit\n'
	tmp += 'fi\n\n'

	tmp += 'mkdir -p $SCRATCH_DIR\n'
	tmp += 'mkdir -p $OUTPUT_DIR\n'
	tmp += 'mkdir -p $INPUT_DIR/ZINDO\n\n'

	tmp += 'cd $INPUT_DIR/MD\n'
	tmp += 'find . -name "*" | cpio -pd $SCRATCH_DIR\n'
	tmp += 'cp $DIR/create_input_zindo* $SCRATCH_DIR\n'
	tmp += 'cd $SCRATCH_DIR\n\n'

	if location_cluster == "lyra" or location_cluster == "adam":
		tmp += '#pgcpp -fast -Minline=levels:10 create_input_zindo.cpp\n'
		tmp += '#mv a.out create_input_zindo\n\n'
	else:
		tmp += 'g++ create_input_zindo.cpp -o create_input_zindo -O2 -lm\n\n'

	tmp += 'for FILE in `find . -name "*.xyz"`\n'
	tmp += 'do\n'
	tmp += '	NAME=`echo $FILE | cut -d "." -f2 | cut -d "/" -f2`\n'
	tmp += '	mkdir -p $NAME\n\n'
	
	tmp += '	i=0\n'
	tmp += "	N_FRAME=`sed -n 1,1p $FILE | awk '{ print $1 }'`\n"
	tmp += '	while [ $i -lt $N_FRAME ]\n'
	tmp += '	do\n'
	tmp += '		mkdir -p $NAME/frame_$i\n'
	tmp += '		mkdir -p output/$NAME/frame_$i\n'
	tmp += '		(( i = $i+1 ))\n'
	tmp += '	done\n\n'

	tmp += '	./create_input_zindo -i $NAME -o output/$NAME -z $ZINDO_DIR\n\n'
	
	tmp += '	i=0\n'
	tmp += '	while [ $i -lt $N_FRAME ]\n'
	tmp += '	do\n'
	tmp += '		tar cfz "$NAME"_frame_"$i".tar.gz $NAME/frame_$i\n'
	tmp += '		mv "$NAME"_frame_"$i".tar.gz $INPUT_DIR/ZINDO\n'
	tmp += '		rm -rf $NAME/frame_$i\n\n'
	
	tmp += '		cd output\n'
	tmp += '		tar cfz output_"$NAME"_frame_"$i"_DIST.tar.gz $NAME/frame_$i\n'
	tmp += '		mv output_"$NAME"_frame_"$i"_DIST.tar.gz $OUTPUT_DIR\n'
	tmp += '		rm -rf output/$NAME/frame_$i\n'
	tmp += '		cd $SCRATCH_DIR\n\n'
	
	tmp += '		(( i = $i+1 ))\n'
	tmp += '	done\n\n'	
	
	tmp += 'done\n\n'

	tmp += 'rm -rf $SCRATCH_DIR/*\n'
	
	try:
		foutput = open('%s/01.file_creation.sh' % (project_name), 'w')
	except:
		print "[ERROR] Could not create %s/01.file_creation.sh. Aborting..." % (project_name)
		sys.exit(1)
	
	foutput.write(tmp)
	foutput.close()
	
def ScriptFileCreationDirect(dir_cluster, input_dir_cluster, output_dir_cluster, zindo_dir_cluster, location_cluster, project_name):
	tmp = ''
	tmp += '#!/bin/bash\n\n'
	tmp += 'DIR="%s"\n' % (dir_cluster)
	tmp += 'INPUT_DIR="%s"\n' % (input_dir_cluster)
	tmp += 'OUTPUT_DIR="%s"\n' % (output_dir_cluster)
	tmp += 'ZINDO_DIR="%s"\n\n' % (zindo_dir_cluster)

	tmp += 'if [[ -d $DIR ]]; then\n'
	tmp += '	cd $DIR\n'
	tmp += 'else\n'
	tmp += '	echo "The folder $DIR does not exist, but it is supposed to be the project directory. Aborting..."\n'
	tmp += '	exit\n'
	tmp += 'fi\n\n'

	tmp += 'mkdir -p $OUTPUT_DIR\n'
	tmp += 'mkdir -p $INPUT_DIR/ZINDO\n\n'

	tmp += 'cp create_input_zindo* $INPUT_DIR/MD\n'
	tmp += 'cd $INPUT_DIR/MD\n\n'

	if location_cluster == "lyra" or location_cluster == "adam":
		tmp += '#module load common pgi\n'
		tmp += '#pgcpp -fast -Minline=levels:10 create_input_zindo.cpp\n'
		tmp += '#mv a.out create_input_zindo\n\n'
	else:
		tmp += 'g++ create_input_zindo.cpp -o create_input_zindo -O2 -lm\n\n'

	tmp += 'for FILE in `find . -name "*.xyz"`\n'
	tmp += 'do\n'
	tmp += '	NAME=`echo $FILE | cut -d "." -f2 | cut -d "/" -f2`\n'
	tmp += '	echo "Currently analyzing file $NAME..."\n'
	tmp += '	mkdir -p $NAME\n\n'
	
	tmp += '	i=0\n'
	tmp += "	N_FRAME=`sed -n 1,1p $FILE | awk '{ print $1 }'`\n"
	tmp += '	while [ $i -lt $N_FRAME ]\n'
	tmp += '	do\n'
	tmp += '		mkdir -p $NAME/frame_$i\n'
	tmp += '		mkdir -p $OUTPUT_DIR/$NAME/frame_$i\n'
	tmp += '		(( i = $i+1 ))\n'
	tmp += '	done\n\n'

	tmp += '	./create_input_zindo -i $NAME -o $OUTPUT_DIR/$NAME -z $ZINDO_DIR\n'

	tmp += '	i=0\n'
	tmp += '	while [ $i -lt $N_FRAME ]\n'
	tmp += '	do\n'
	tmp += '		cd $INPUT_DIR/MD\n'
	tmp += '		tar cfz "$NAME"_frame_"$i".tar.gz $NAME/frame_$i\n'
	tmp += '		mv "$NAME"_frame_"$i".tar.gz $INPUT_DIR/ZINDO\n'
	tmp += '		rm -rf $NAME/frame_$i\n\n'
	
	tmp += '		cd $OUTPUT_DIR\n'
	tmp += '		tar cfz output_"$NAME"_frame_"$i"_DIST.tar.gz $NAME/frame_$i\n'
	tmp += '		rm -rf $OUTPUT_DIR/$NAME/frame_$i\n\n'
	
	tmp += '		(( i = $i+1 ))\n'
	tmp += '	done\n\n'
	tmp += '	rm -rf $OUTPUT_DIR/$NAME\n\n'	
	tmp += 'done\n\n'
	tmp += 'rm $INPUT_DIR/MD/create_input_zindo*\n'
	
	try:
		foutput = open('%s/01.file_creation_direct.sh' % (project_name), 'w')
	except:
		print "[ERROR] Could not create %s/01.file_creation_direct.sh. Aborting..." % (project_name)
		sys.exit(1)
	
	foutput.write(tmp)
	foutput.close()

def ScriptFileCreationPBS(project_name, username_cluster, dir_cluster, location_cluster):
	tmp = ''
	if location_cluster == "joe":
		tmp += '#!/bin/csh\n'
		tmp += '#PBS -q long\n' 
		tmp += '#PBS -l nodes=1:p4:ppn=1\n'
		tmp += '#PBS -A %s\n' % (username_cluster)
		tmp += '#PBS -M %s@averell.umh.ac.be\n' % (username_cluster)
		tmp += '#PBS -m bae\n'
		tmp += '#PBS -V\n\n'
	
	elif location_cluster == "lyra" or location_cluster == "adam":
		tmp += '#!/bin/bash\n\n'
		tmp += '#$ -j y\n'
		tmp += '#$ -cwd\n'
		tmp += '#$ -l vf=2G\n'
		tmp += '#$ -l h_cpu=600:00:00\n'
		tmp += '#$ -N %s\n' % (project_name)
		tmp += '#$ -m bea\n'
		tmp += '#$ -M nicolas.g.martinelli@gmail.com\n\n'
 
		tmp += 'module load common pgi\n\n'
	
	else:
		print '[ERROR] Bad cluster location. Aborting...'
		sys.exit(1)

	tmp += 'cd %s\n' % (dir_cluster)
	tmp += 'chmod +x 01.file_creation.sh\n'
	tmp += './01.file_creation.sh\n'

	try:
		foutput = open('%s/01.file_creation.pbs' % (project_name), 'w')
	except:
		print "[ERROR] Could not create %s/01.file_creation.pbs. Aborting..." % (project_name)
		sys.exit(1)
	
	foutput.write(tmp)
	foutput.close()

def ScriptZINDOLaunch(dir_cluster, input_dir_cluster, output_dir_cluster, scratch_dir_cluster, project_name, location_cluster):
	
	tmp = ''
	tmp += '#!/bin/bash\n\n'
	tmp += 'DIR="%s"\n' % (dir_cluster)
	tmp += 'INPUT_DIR="%s"\n' % (input_dir_cluster)
	tmp += 'OUTPUT_DIR="%s"\n' % (output_dir_cluster)
	if location_cluster == "lyra" or location_cluster == "adam":
		tmp += 'SCRATCH_DIR="\\%s"\n\n' % (scratch_dir_cluster)
	else:
		tmp += 'SCRATCH_DIR="%s"\n\n' % (scratch_dir_cluster)
	tmp += 'N_PBS=6\n\n'

	tmp += 'MakePBS(){\n'

	if location_cluster == "joe":
		tmp += '	echo "#!/bin/csh" 				> $DIR/zindo_$1.pbs\n'
		tmp += '	echo "#PBS -q long"				>> $DIR/zindo_$1.pbs\n'
		tmp += '	echo "#PBS -l nodes=1:p4:ppn=1" 		>> $DIR/zindo_$1.pbs\n'
		tmp += '	echo "#PBS -A `whoami`" 			>> $DIR/zindo_$1.pbs\n'
		tmp += '	echo "#PBS -M `whoami`@averell.umh.ac.be" 		>> $DIR/zindo_$1.pbs\n'
		tmp += '	echo "#PBS -m bae" 				>> $DIR/zindo_$1.pbs\n'
		tmp += '	echo "#PBS -V" 					>> $DIR/zindo_$1.pbs\n'
		tmp += '	echo " "					>> $DIR/zindo_$1.pbs\n'

	elif location_cluster == "lyra" or location_cluster == "adam":
		tmp += '	echo "#!/bin/bash"			>> $DIR/zindo_$1.pbs\n'
		tmp += '	echo "#$ -j y"					>> $DIR/zindo_$1.pbs\n'
		tmp += '	echo "#$ -cwd"					>> $DIR/zindo_$1.pbs\n'
		tmp += '	echo "#$ -l vf=2G"				>> $DIR/zindo_$1.pbs\n'
		tmp += '	echo "#$ -l h_cpu=600:00:00"		>> $DIR/zindo_$1.pbs\n'
		tmp += '	echo "#$ -N %s"		>> $DIR/zindo_$1.pbs\n' % (project_name)
		tmp += '	echo "#$ -m bea"		>> $DIR/zindo_$1.pbs\n'
		tmp += '	echo "#$ -M nicolas.g.martinelli@gmail.com"		>> $DIR/zindo_$1.pbs\n'
 		tmp += '	echo " "					>> $DIR/zindo_$1.pbs\n'
		tmp += '	echo "module load common"			>> $DIR/zindo_$1.pbs\n'
		tmp += '	echo " "					>> $DIR/zindo_$1.pbs\n'
	
	else:
		print '[ERROR] Bad cluster location. Aborting...'
		sys.exit(1)
		
	tmp += '	echo "cd $DIR"					>> $DIR/zindo_$1.pbs\n'
	tmp += '	echo "chmod +x zindo_$1.run"			>> $DIR/zindo_$1.pbs\n'
	tmp += '	echo "./zindo_$1.run"				>> $DIR/zindo_$1.pbs\n'
	tmp += '}\n\n'

#	tmp += 'MakeRUN(){\n'
#	tmp += '	echo "#!/bin/bash" 				> $DIR/zindo_$1.run\n'
#	tmp += '	echo " "					>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "mkdir -p $SCRATCH_DIR"			>> $DIR/zindo_$1.run\n'
#	tmp += '	echo " "					>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "for x in \`cat $DIR/zindo_$1.dir\`"	>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "do"					>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	cd $INPUT_DIR/ZINDO"			>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	cp --parents -R -f \$x $SCRATCH_DIR"	>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	cd $SCRATCH_DIR/\$x"			>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	for CMD in \`find . -name \'*.cmd\'\`"	>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	do"					>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "		chmod +x \$CMD"			>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "		./\$CMD"			>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	done"					>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	rm -rf $SCRATCH_DIR/\$x"		>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "done"					>> $DIR/zindo_$1.run\n\n'
#	tmp += '}\n\n'

#	tmp += 'MakeRUN(){\n'
#	tmp += '	echo "#!/bin/bash"								> $DIR/zindo_$1.run\n'
#	tmp += '	echo " "										>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "mkdir -p $SCRATCH_DIR"					>> $DIR/zindo_$1.run\n'
#	tmp += '	echo " "										>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "for TAR in \`cat $DIR/zindo_$1.dir\`"		>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "do"										>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	cd $INPUT_DIR/ZINDO"					>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	cp \$TAR $SCRATCH_DIR"					>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	cd $SCRATCH_DIR"						>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	x=\`echo \$TAR | cut -d \"/\" -f2 | cut -d \".\" -f1\`"	>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	mkdir \$x ; mv \$TAR \$x ; cd \$x"		>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	tar xfz \$TAR"							>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	for FILE in \`find . -name \'*.*\'\`"		>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	do"										>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "		mv \$FILE ."					>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	done"									>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	for CMD in \`find . -name \'*.cmd\'\`"	>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	do"										>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "		chmod +x \$CMD"					>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "		./\$CMD"						>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	done"									>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "	rm -rf $SCRATCH_DIR/\$x"				>> $DIR/zindo_$1.run\n'
#	tmp += '	echo "done"										>> $DIR/zindo_$1.run\n'

	tmp += 'MakeRUN(){\n'
	tmp += '	echo "#!/bin/bash" > $DIR/zindo_$1.run\n\n'
	tmp += '	echo " " >> $DIR/zindo_$1.run\n'
	
	tmp += '	echo "OUTPUT_DIR="%s"" >> $DIR/zindo_$1.run\n\n' % (output_dir_cluster)
	tmp += '	echo " " >> $DIR/zindo_$1.run\n'
	
	tmp += '	echo "mkdir -p $SCRATCH_DIR" >> $DIR/zindo_$1.run\n\n'
	tmp += '	echo " " >> $DIR/zindo_$1.run\n'
	
	tmp += '	echo "for TAR in \`cat $DIR/zindo_$1.dir\`" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "do" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "	cd $INPUT_DIR/ZINDO" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "	cp \$TAR $SCRATCH_DIR" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "	cd $SCRATCH_DIR" >> $DIR/zindo_$1.run\n\n'
	tmp += '	echo " " >> $DIR/zindo_$1.run\n'
	
	tmp += '	echo "	PROJECT=\`echo \$TAR | awk -F \'/\' \'{print \$2}\' | awk -F \'_frame_\' \'{print \$1}\'\`" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "	FRAME=\`echo \$TAR | awk -F \'_frame_\' \'{print \$2}\' | cut -d \'.\' -f1\`" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "	mkdir -p \$PROJECT/frame_\$FRAME" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "	mv \$TAR \$PROJECT/frame_\$FRAME" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "	cd \$PROJECT/frame_\$FRAME" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "	tar xfz \$TAR" >> $DIR/zindo_$1.run\n\n'
	tmp += '	echo " " >> $DIR/zindo_$1.run\n'
	
	tmp += '	echo "	for FILE in \`find . -name \'*.*\'\`" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "	do" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "		mv \$FILE ." >> $DIR/zindo_$1.run\n'
	tmp += '	echo "	done" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "	for CMD in \`find . -name \'*.cmd\'\`" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "	do" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "		chmod +x \$CMD" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "		./\$CMD" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "	done" >> $DIR/zindo_$1.run\n\n'
	tmp += '	echo " " >> $DIR/zindo_$1.run\n'
	
	tmp += '	echo "	cd $SCRATCH_DIR" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "	tar cfz output_\\"\$PROJECT\\"_frame_\\"\$FRAME\\"_J.tar.gz \$PROJECT/frame_\$FRAME/dimer_*_*.out" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "	mv output_\\"\$PROJECT\\"_frame_\\"\$FRAME\\"_J.tar.gz \$OUTPUT_DIR" >> $DIR/zindo_$1.run\n\n'
	tmp += '	echo " " >> $DIR/zindo_$1.run\n'
	
	tmp += '	echo "	rm -rf \$PROJECT/frame_\$FRAME" >> $DIR/zindo_$1.run\n'
	tmp += '	echo "done" >> $DIR/zindo_$1.run\n'
	tmp += '}\n\n'
	
	tmp += 'i=1\n'
	tmp += 'j=1\n'
	tmp += 'k=0\n\n'

	tmp += 'cd $INPUT_DIR/ZINDO\n\n'

	tmp += 'find . -name "*frame*" > directories.tmp\n'
	tmp += 'N_DIR=`wc -l directories.tmp | awk \'{print $1}\'`\n'
	tmp += 'N_STEP=$(($N_DIR/$N_PBS + 1))\n'

	tmp += 'while [ $i -le $N_PBS ]\n'
	tmp += 'do\n\n'

	tmp += '	k=$(($j+$N_STEP))\n'
	tmp += '	MakePBS $i\n'
	tmp += '	MakeRUN $i\n'
	tmp += '	sed -n "$j","$k"p directories.tmp > $DIR/zindo_$i.dir\n\n'

	tmp += '	j=$(($k+1))\n'
	tmp += '	i=$(($i+1))\n\n'

	tmp += 'done\n\n'
	
	tmp += 'cd $DIR\n'
	tmp += 'for PBS in `ls zindo_*.pbs`\n'
	tmp += 'do\n'
	tmp += '	qsub $PBS\n'
	tmp += 'done\n'

	
	try:
		foutput = open('%s/02.launch_zindo.sh' % (project_name), 'w')
	except:
		print "[ERROR] Could not create %s/02.launch_zindo.sh. Aborting..." % (project_name)
		sys.exit(1)
	
	foutput.write(tmp)
	foutput.close()
