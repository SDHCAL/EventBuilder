 if [ "$#" -eq  "0" ]
   then
 echo "all? streamout? Trivent? Analysis? ... "
else 

  if [[ $1 == all ]]
  then
  for n in DHCAL_*_I0_0.slcio
  do
  m=$( echo "$n" | cut -c 7-12 ) 
	touch DHCAL_"$m".txt
	for l in DHCAL_"$m"_I0_*.slcio
	do
		echo $(pwd)/$l>>DHCAL_"$m".txt
	done
	sed -re "12r DHCAL_"$m".txt" "./xml/Streamout.xml">DHCAL_"$m".xml
	rm DHCAL_"$m".txt
	sed -i "s|ModifieMoi|$(pwd)/DHCAL_Streamout_"$m"_I0.slcio|" DHCAL_"$m".xml
	cp $(pwd)/xml/Trivent.xml $(pwd)
	sed -i "s|FILE|$(pwd)/DHCAL_Streamout_"$m"_I0.slcio|" Trivent.xml
	sed -i "s|OUTPUT|$(pwd)/DHCAL_Trivent_"$m"_I0.slcio|" Trivent.xml
	sed -i "s|NOISES|$(pwd)/DHCAL_Noise_"$m"_I0.slcio|" Trivent.xml
	cp $(pwd)/xml/Analysis.xml $(pwd)/Analysis_"$m".xml
	sed -i "s|FILE|$(pwd)/DHCAL_Trivent_"$m"_I0.slcio|" Analysis_"$m".xml
	Marlin DHCAL_"$m".xml
	rm DHCAL_"$m".xml
	Marlin Trivent.xml
	rm Trivent.xml
	Marlin Analysis_"$m".xml
	rm Analysis_"$m".xml
	done
	elif [[ $1 == streamout ]]
	then
	for n in DHCAL_*_I0_0.slcio
	do
	m=$( echo "$n" | cut -c 7-12 ) 
	touch DHCAL_"$m".txt
	for l in DHCAL_"$m"_I0_*.slcio
	do
		echo $(pwd)/$l>>DHCAL_"$m".txt
	done
	sed -re "12r DHCAL_"$m".txt" "./xml/Streamout.xml">DHCAL_"$m".xml
	rm DHCAL_"$m".txt
        if [ "$#" -eq  "1" ] 
        then sed -i "s|NUMBER|0|" DHCAL_"$m".xml
        else sed -i "s|NUMBER|$2|" DHCAL_"$m".xml
        fi
	sed -i "s|OUTPUT|$(pwd)/DHCAL_Streamout_"$m"_I0.slcio|" DHCAL_"$m".xml
	Marlin DHCAL_"$m".xml
	rm DHCAL_"$m".xml
	done
	elif [[ $1 == analysis ]]
	then
	for n in DHCAL_Trivent_*_I0.slcio
	do
	m=$( echo "$n" | cut -c 15-20 ) 
	cp $(pwd)/xml/Analysis.xml $(pwd)/Analysis_"$m".xml
        if [ "$#" -eq  "1" ] 
        then sed -i "s|NUMBER|0|" Analysis_"$m".xml
        else sed -i "s|NUMBER|$2|" Analysis_"$m".xml
        fi
	sed -i "s|FILE|$(pwd)/DHCAL_Trivent_"$m"_I0.slcio|" Analysis_"$m".xml
	Marlin Analysis_"$m".xml
	rm Analysis_"$m".xml
        now=$(date +"%m_%d_%Y")
        #mv Results_"$m".txt $(pwd)/Results_"$m"_"$now".txt
        #mv Results_Analysis_"$m".root $(pwd)/Results_Analysis_"$m"_"$now".root
	done
        elif [[ $1 == analysisSplit ]]
	then
	for n in DHCAL_Trivent_*_I0_Split.slcio
	do
	m=$( echo "$n" | cut -c 15-20 ) 
	cp $(pwd)/xml/AnalysisSplit.xml $(pwd)/Analysis_"$m".xml
	sed -i "s|FILE|$(pwd)/DHCAL_Trivent_"$m"_I0_Split.slcio|" Analysis_"$m".xml
	Marlin Analysis_"$m".xml >>Result_"$m"_Split.txt
	rm Analysis_"$m".xml
        now=$(date +"%m_%d_%Y")
        #mv Results_"$m".txt $(pwd)/Split/Results_"$m"_"$now".txt
        #mv Results_Analysis_"$m".root $(pwd)/Split/Results_Analysis_"$m"_"$now".root
	done
        mv Results.txt $(pwd)/Split/Results_Split_"$now".txt
        elif [[ $1 == noise ]]
	then
	for n in DHCAL_Noise_*_I0.slcio
	do
	m=$( echo "$n" | cut -c 13-18 ) 
	cp $(pwd)/xml/Noise.xml $(pwd)/Noise_"$m".xml
        if [ "$#" -eq  "1" ] 
        then sed -i "s|NUMBER|0|" Noise_"$m".xml
        else sed -i "s|NUMBER|$2|" Noise_"$m".xml
        fi
	sed -i "s|FILE|$(pwd)/DHCAL_Noise_"$m"_I0.slcio|" Noise_"$m".xml
	Marlin Noise_"$m".xml
	rm Noise_"$m".xml
	done
	elif [[ $1 == trivent ]]
	then
	for n in DHCAL_Streamout_*_I0.slcio
	do
	cp $(pwd)/xml/Trivent.xml $(pwd)
	m=$( echo "$n" | cut -c 17-22 )
        if [ "$#" -eq  "1" ] 
        then sed -i "s|NUMBER|0|" Trivent.xml
        else sed -i "s|NUMBER|$2|" Trivent.xml
        fi
	sed -i "s|FILE|$(pwd)/DHCAL_Streamout_"$m"_I0.slcio|" Trivent.xml
	sed -i "s|OUTPUT|$(pwd)/DHCAL_Trivent_"$m"_I0.slcio|" Trivent.xml
	sed -i "s|NOISES|$(pwd)/DHCAL_Noise_"$m"_I0.slcio|" Trivent.xml
        now=$(date +"%m_%d_%Y")
        #mv Results_"m".root $(pwd)/Results_Trivent_"m"_"$now".root
	Marlin Trivent.xml
	rm Trivent.xml
	done
        elif [[ $1 == triventnonoise ]]
	then
	for n in DHCAL_Streamout_*_I0.slcio
	do
	cp $(pwd)/xml/Trivent.xml $(pwd)
	m=$( echo "$n" | cut -c 17-22 ) 
        if [ "$#" -eq  "1" ] 
        then sed -i "s|NUMBER|0|" Trivent.xml
        else sed -i "s|NUMBER|$2|" Trivent.xml
        fi
	sed -i "s|FILE|$(pwd)/DHCAL_Streamout_"$m"_I0.slcio|" Trivent.xml
	sed -i "s|OUTPUT|$(pwd)/DHCAL_Trivent_"$m"_I0.slcio|" Trivent.xml
	sed -i "s|NOISES| |" Trivent.xml
        now=$(date +"%m_%d_%Y")
        #mv Results_"m".root $(pwd)/Results_Trivent_"m"_"$now".root
	Marlin Trivent.xml
	rm Trivent.xml
	done
        elif [[ $1 == triventSplit ]]
	then
	for n in DHCAL_Streamout_*_I0.slcio
	do
	cp $(pwd)/xml/TriventSplit.xml $(pwd)
	m=$( echo "$n" | cut -c 17-22 ) 
	sed -i "s|FILE|$(pwd)/DHCAL_Streamout_"$m"_I0.slcio|" TriventSplit.xml
	sed -i "s|OUTPUT|$(pwd)/DHCAL_Trivent_"$m"_I0_Split.slcio|" TriventSplit.xml
	sed -i "s|NOISES|$(pwd)/DHCAL_Noise_"$m"_I0_Split.slcio|" TriventSplit.xml
	Marlin TriventSplit.xml
        #mv Results_Trivent_*.root $(pwd)/Split
	rm TriventSplit.xml
	done
	fi
fi

