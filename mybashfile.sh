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
	sed -re "12r DHCAL_"$m".txt" "./xml/streamout.xml">DHCAL_"$m".xml
	rm DHCAL_"$m".txt
	sed -i "s|ModifieMoi|$(pwd)/DHCAL_Streamout_"$m"_I0.slcio|" DHCAL_"$m".xml
	cp $(pwd)/xml/Trivent.xml $(pwd)
	sed -i "s|Salut|$(pwd)/DHCAL_Streamout_"$m"_I0.slcio|" Trivent.xml
	sed -i "s|ModifieMoi|$(pwd)/DHCAL_Trivent_"$m"_I0.slcio|" Trivent.xml
	sed -i "s|Test|$(pwd)/DHCAL_Noise_"$m"_I0.slcio|" Trivent.xml
	cp $(pwd)/xml/Analysis.xml $(pwd)/Analysis_"$m".xml
	sed -i "s|Salut|$(pwd)/DHCAL_Trivent_"$m"_I0.slcio|" Analysis_"$m".xml
	sed -i "s|yuii|"$m"|" Analysis_"$m".xml
	Marlin DHCAL_"$m".xml
	rm DHCAL_"$m".xml
	Marlin Trivent.xml
	rm Trivent.xml
	Marlin Analysis_"$m".xml >>Result_"$m".txt
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
	sed -re "12r DHCAL_"$m".txt" "./xml/streamout.xml">DHCAL_"$m".xml
	rm DHCAL_"$m".txt
	sed -i "s|ModifieMoi|$(pwd)/DHCAL_Streamout_"$m"_I0.slcio|" DHCAL_"$m".xml
	Marlin DHCAL_"$m".xml
	rm DHCAL_"$m".xml
	done
	elif [[ $1 == analysis ]]
	then
	for n in DHCAL_Trivent_*_I0.slcio
	do
	m=$( echo "$n" | cut -c 15-20 ) 
	cp $(pwd)/xml/Analysis.xml $(pwd)/Analysis_"$m".xml
	sed -i "s|Salut|$(pwd)/DHCAL_Trivent_"$m"_I0.slcio|" Analysis_"$m".xml
	sed -i "s|yuii|"$m"|" Analysis_"$m".xml
	Marlin Analysis_"$m".xml >>Result_"$m".txt
	rm Analysis_"$m".xml
	done
        elif [[ $1 == noise ]]
	then
	for n in DHCAL_Noise_*_I0.slcio
	do
	m=$( echo "$n" | cut -c 13-18 ) 
	cp $(pwd)/xml/Noise.xml $(pwd)/Noise_"$m".xml
	sed -i "s|Salut|$(pwd)/DHCAL_Noise_"$m"_I0.slcio|" Noise_"$m".xml
	sed -i "s|yuii|"$m"|" Noise_"$m".xml
	Marlin Noise_"$m".xml
	rm Noise_"$m".xml
	done
	elif [[ $1 == trivent ]]
	then
	for n in DHCAL_Streamout_*_I0.slcio
	do
	cp $(pwd)/xml/Trivent.xml $(pwd)
	m=$( echo "$n" | cut -c 17-22 ) 
	sed -i "s|Salut|$(pwd)/DHCAL_Streamout_"$m"_I0.slcio|" Trivent.xml
	sed -i "s|ModifieMoi|$(pwd)/DHCAL_Trivent_"$m"_I0.slcio|" Trivent.xml
	sed -i "s|Test|$(pwd)/DHCAL_Noise_"$m"_I0.slcio|" Trivent.xml
	Marlin Trivent.xml
	rm Trivent.xml
	done
	fi
        if [[ $1 == pdf ]]
 	then
        $(pwd)/foo $2 $3
        fi
fi

