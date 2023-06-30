#*##########################################################################*#
##################        Blip Analysis Functions           ##################
#*##########################################################################*#

def checkBlipSDev(df_to_test): #Calculate standard deviation of blip charges in each event in the pandas dataframe
    df_to_test_grouped = df_to_test.groupby(["event","run","timestamp"])
    goodEventNum = 0
    numEvents = 0
    for group_name, df_test_group in df_to_test_grouped:
        eventQs = []
        for q in df_test_group["blip_cluster_charge"]:
            eventQs.append(q)
        if len(eventQs) > 1:
            qStdev = np.std(eventQs)
            #print "Sum blip q: ", np.sum(eventQs)," \tMean blip q: ",np.mean(eventQs), "\tStdev q's: ",qStdev
            if qStdev > 1.01:
                goodEventNum += 1

        numEvents += 1
    print "Number events out of %s with variable blip_q: "%(numEvents), goodEventNum
    print "Percent good events: ", 100*goodEventNum/numEvents

def countBlipMultiplicity(df_to_test,complete_df): #counts blip multiplicity around a max point
    numMaxPts = 6 #max number of clusters to sum
    multiplicityInRadius = []
    summedBlipEs = []

    df_grouped = df_to_test.groupby(["event","run","timestamp"])
    #cd_grouped = completeDetDF(["event","run","timestamp"])
    print "Number of events with energies: ", df_grouped["blip_energy"].first().count() #counting number of unique events
    counter = 0

    for group_name, df_group in df_grouped:
        usedEnergies = []
        blip_Es_ievent = df_group["blip_energy"].agg(list)
        maxesUsed = 0 #Counts how many of the higher energy blips in an event were used. Max currently 6
        for j in range(len(df_group)):
            if len(blip_Es_ievent) > 0:
                curr_Max_E = max(blip_Es_ievent)
                if curr_Max_E < 0.45:    #Only want E > 0.500 MeV
                    break
                maxesUsed+=1
                multiplicity = 1
                maxIndex = blip_Es_ievent.index(curr_Max_E)
                usedEnergies.append(curr_Max_E)
                blip_Es_ievent.remove(curr_Max_E)
                maxMember = df_group.iloc[maxIndex]
                blip_x,blip_y,blip_z = maxMember["blip_x"],maxMember["blip_y"],maxMember["blip_z"]
                nRun,nSubRun,nEvent = maxMember["run"],maxMember["timestamp"],maxMember["event"]

                nearby_blips = complete_df.query('(run == %s) and (timestamp == %s) and(event == %s) and (blip_x > %s and blip_x < %s) and (blip_y > %s and blip_y < %s) and (blip_z>%s and blip_z< %s)'%(nRun, nSubRun, nEvent,blip_x-20.0,blip_x+20.0,blip_y-20.0,blip_y+20.0,blip_z-20.0,blip_z+20.0))
                completeNearbyEs = nearby_blips["blip_energy"].agg(list)

                nearbyESums = curr_Max_E
                for nearE in completeNearbyEs:
                    if nearE not in usedEnergies:
                        usedEnergies.append(nearE)
                        nearbyESums += nearE
                        multiplicity += 1
                        if nearE in blip_Es_ievent:
                            blip_Es_ievent.remove(nearE)
                multiplicityInRadius.append(multiplicity)
                summedBlipEs.append(nearbyESums)
            elif maxesUsed == numMaxPts:
                break

        counter+=1
        if counter > 3000:
                break
    print "Number of energy clusters in this DF: ", len(summedBlipEs), " out of "," energies total"
    return summedBlipEs,multiplicityInRadius

def DistTweenTwoPoints(memberA,memberB): #calculates the real 3D distance between two blips or events in the detector
    xA = memberA["blip_x"]
    yA = memberA["blip_y"]
    zA = memberA["blip_z"]
    xB = memberB["blip_x"]
    yB = memberB["blip_y"]
    zB = memberB["blip_z"]
    deltaxSqd = ((xA-xB)**2)
    deltaySqd = ((yA-yB)**2)
    deltazSqd= ((zA-zB)**2)
    dist_from_A = np.sqrt(deltaxSqd+deltaySqd+deltazSqd)
    #print "x1,x2: (%s,%s)\ty1,y2: (%s,%s)\tz1,z2: (%s,%s)\t"%(xA,xB,yA,yB,zA,zB)
    #print "Distance: ",dist_from_A
    return dist_from_A

def DistTweenMaxNPoint(xA,yA,zA,memberB):
    xB = memberB["blip_x"]
    yB = memberB["blip_y"]
    zB = memberB["blip_z"]
    deltaxSqd = ((xA-xB)**2)
    deltaySqd = ((yA-yB)**2)
    deltazSqd= ((zA-zB)**2)
    dist_from_A = np.sqrt(deltaxSqd+deltaySqd+deltazSqd)
    #print "x1,x2: (%s,%s)\ty1,y2: (%s,%s)\tz1,z2: (%s,%s)\t"%(xA,xB,yA,yB,zA,zB)
    #print "Distance: ",dist_from_A
    return dist_from_A

def distancesFromPt(Loc,thisGroup):#distance between one blip and all other blips in event
    distances = []
    deltaxSqds = []
    deltaySqds = []
    deltazSqds = []
    xMaxblipQ = Loc["blip_x"]
    yMaxblipQ = Loc["blip_y"]
    zMaxblipQ = Loc["blip_z"]
    #print "x: ",xMaxblipQ," y: ",yMaxblipQ," z: ",zMaxblipQ," "
    #print xMaxblipQ," ",yMaxblipQ," ",zMaxblipQ
    for x in thisGroup["blip_x"]:
        deltaxSqds.append((xMaxblipQ-x)**2)
    for y in thisGroup["blip_y"]:
        deltaySqds.append((yMaxblipQ-y)**2)
    for z in thisGroup["blip_z"]:
        deltazSqds.append((zMaxblipQ-z)**2)
    for i in range(len(deltaxSqds)):
        dist_from_m = np.sqrt(deltaxSqds[i]+deltaySqds[i]+deltazSqds[i])
        #print "distfromm: ",dist_from_m
        if dist_from_m > 0.05: #If not the max point...
            distances.append(np.sqrt(deltaxSqds[i]+deltaySqds[i]+deltazSqds[i]))
        else:
            distances.append(50000)
    return distances

def distancesFromPtSlick(ptLoc,thisGroup,NumRows): #Shorter, but actually slower by factor of 10)
    f_distances = []
    for k in range(NumRows):
        rowk = thisGroup.iloc[k]
        dist_from_pt = np.sqrt((((ptLoc["blip_x"])-(rowk["blip_x"]))**2)+(((ptLoc["blip_y"])-(rowk["blip_y"]))**2)+(((ptLoc["blip_z"])-(rowk["blip_z"]))**2))
        if dist_from_pt > 0.05: #If not the max point...
            f_distances.append(dist_from_pt)
        else:
            f_distances.append(50000)
    return f_distances

def fractionAboveEnergy(thisDF,E_lower_bound):  #Calculate the percentage of blips above a certain minimum energy
    #print "working location"
    all_energies = thisDF["blip_energy"]
    thisDF_E_cut = thisDF.query('energy>%s'%(E_lower_bound))
    cut_energies = thisDF_E_cut["blip_energy"]
    print "Fraction: ", (float(len(cut_energies))/float(len(all_energies))),"\n From ", len(cut_energies), " out of total ", len(all_energies)

def fractionEnergiesWithNoNeutrinos(thisDF):  #Find the percent of events with energies with no neutrinos in them
    #Fraction for blip branch... should do for event branch
    cNoNu = 'neutrinos == 0 and neutrinoshowers==0 and cosmic_trk_50==0'
    #print "working location"
    all_energies = thisDF["blip_energy"]
    thisDF_Nu_cut = thisDF.query(cNoNu)
    cut_energies = thisDF_Nu_cut["blip_energy"]
    print "Fraction: ", (float(len(cut_energies))/float(len(all_energies))),"\n From ", len(cut_energies), " out of total ", len(all_energies)

def fractionWithNoNeutrinos(thisDF): #Find the percent of events with no neutrinos in them
    #Fraction for blip branch... should do for event branch
    cNoNu = 'neutrinos == 0 and neutrinoshowers==0 and cosmic_trk_50==0'
    #print "working location"
    all_energies = thisDF["blip_energy"]
    thisDF_Nu_cut = thisDF.query(cNoNu)
    cut_energies = thisDF_Nu_cut["blip_energy"]
    print "Fraction: ", (float(len(cut_energies))/float(len(all_energies))),"\n From ", len(cut_energies), " out of total ", len(all_energies)

def fractionEVTsWithNoNeutrinos(thisDF_EVT):
    evt_times = thisDF_EVT["evttime"]
    thisDF_Nu_cut = thisDF_EVT.query(cNoNu)
    cut_evt_ts = thisDF_Nu_cut["evttime"]
    print "Fraction events with no neutrinos: ", (float(len(cut_evt_ts ))/float(len(evt_times))),"\n Or ", len(cut_evt_ts), " out of total ", len(evt_times)

def groupDataEventsAndAddMaxblip(thisDF):
    summedCloseClusters = [] #for max blip_Q in an event + nearest neighbor
    maxblips = []
    df_grouped = thisDF.groupby(["event","run","timestamp"])
    print "Number of max blip points: ", df_grouped["blip_cluster_charge"].max().count() #counting number of unique events
    #Here, we're cycling through all the stuff in one event -- since we're on blip branch, this is a collection of blip points
    for group_name, df_group in df_grouped:
        #print df_group[group_name]
        #print df_group
        #print "\nDf_group EV min: ",df_group["blip_cluster_charge"].min()

        eventMaxq = df_group["blip_cluster_charge"].max()
        #print "EventMaxq: ",eventMaxq
        if eventMaxq > -0.01:
            maxLocation = df_group.loc[df_group["blip_cluster_charge"].idxmax()]
        else:
            continue

        blip_qs_ievent = df_group["blip_cluster_charge"].agg(list)
        distancesfromMax = distancesFromPt(maxLocation,df_group)
        #print "d's: ",distancesfromMax

        if len(distancesfromMax) != len(blip_qs_ievent):
            print "Length mismatch... in distances and q's"
        else:
            maxblips.append(eventMaxq)
            qindex = distancesfromMax.index(min(distancesfromMax))
            summedCloseClusters.append(eventMaxq+blip_qs_ievent[qindex])
            #print "Sum q's: ",summedCloseClusters
    return maxblips,summedCloseClusters


def groupDataEventsAndAddMaxblipConditional(thisDF):
    summedCloseClusters = [] #for max blip_Q in an event + nearest neighbor
    maxblips = []
    df_grouped = thisDF.groupby(["event","run","timestamp"])
    print "Number of max blip points: ", df_grouped["blip_cluster_charge"].max().count() #counting number of unique events
    #Here, we're cycling through all the stuff in one event -- since we're on blip branch, this is a collection of blip points
    for group_name, df_group in df_grouped:
        #print df_group[group_name]
        #print df_group
        #print "\nDf_group EV min: ",df_group["blip_cluster_charge"].min()

        eventMaxq = df_group["blip_cluster_charge"].max()
        #print "EventMaxq: ",eventMaxq
        if eventMaxq > -0.01:
            maxLocation = df_group.loc[df_group["blip_cluster_charge"].idxmax()]
        else:
            continue

        blip_qs_ievent = df_group["blip_cluster_charge"].agg(list)
        distancesfromMax = distancesFromPt(maxLocation,df_group)
        print "d's: ",distancesfromMax

        if len(distancesfromMax) != len(blip_qs_ievent):
            print "Length mismatch... in distances and q's"
        else:
            maxblips.append(eventMaxq)
            qindex = distancesfromMax.index(min(distancesfromMax))
            summedCloseClusters.append(eventMaxq+blip_qs_ievent[qindex])
            #print "Sum q's: ",summedCloseClusters
    return maxblips,summedCloseClusters

def returnMaxBlipEs(thisDF): #Find the highest blip energy in an event, or eventdisplay
    maxblips = []
    df_grouped = thisDF.groupby(["event","run","timestamp"])
    #Here, we're cycling through all the stuff in one event -- since we're on blip branch, this is a collection of blip points
    for group_name, df_group in df_grouped:
        eventMaxE = df_group["blip_energy"].max()
        maxblips.append(eventMaxE)
    return maxblips

def plotMaxBlipLocs(thisDF):  #plot in two dimensions the location of the highest energy blip deposite in a data event
    #print "working location"
    thisDF_cut = thisDF.query('energy>2.3' )
    fig = plt.figure(figsize=(12,7))
    plt.hist2d(thisDF_cut['blip_z'].values,thisDF_cut['blip_y'].values,bins=(100,100), range=[[-100,1100], [-120,120]],vmin=0,label='OFF Beam Data')#,norm=LogNorm())
    plt.colorbar()
    plt.title("blip, Events with E > 2.3 MeV")
    plt.xlabel('Z')
    plt.ylabel('Y')

plt.tight_layout()

def returnMaxBlipQs(thisDF): #return an array of the max blip charges in the dataframe
    maxblips = []
    df_grouped = thisDF.groupby(["event","run","timestamp"])
    print "Number of max blip points: ", df_grouped["blip_cluster_charge"].max().count() #counting number of unique events
    #Here, we're cycling through all the stuff in one event -- since we're on blip branch, this is a collection of blip points
    for group_name, df_group in df_grouped:
        eventMaxq = df_group["blip_cluster_charge"].max()
        maxblips.append(eventMaxq)
    return maxblips

def returnMaxBlipsPlusNeighbors(thisDF,completeDetDF,summingRadius):  #Return the max blip energies + nearby neighbors for total energy in a region
    numMaxPts = 4 #max number of clusters to sum

    summedBlipEs = [] #for First blip_Q in an event + nearest neighbor + next nearest neighbor
    multiplicityInRadius = []

    df_grouped = thisDF.groupby(["event","run","timestamp"])
    #cd_grouped = completeDetDF(["event","run","timestamp"])
    print "Number of events with energies: ", df_grouped["blip_energy"].first().count() #counting number of unique events
    counter = 0

    for group_name, df_group in df_grouped:
        usedEnergies = []
        blip_Es_ievent = df_group["blip_energy"].agg(list)
        maxesUsed = 0 #Counts how many of the higher energy blips in an event were used. Max currently 3
        for j in range(len(df_group)):
            if len(blip_Es_ievent) > 0 and maxesUsed < numMaxPts:
                curr_Max_E = max(blip_Es_ievent)
                if curr_Max_E < 0.50:    #Only want E > 0.500 MeV
                    break
                maxesUsed+=1
                multiplicity = 1
                maxIndex = blip_Es_ievent.index(curr_Max_E)
                usedEnergies.append(curr_Max_E)
                blip_Es_ievent.remove(curr_Max_E)
                maxMember = df_group.iloc[maxIndex]
                blip_x,blip_y,blip_z = maxMember["blip_x"],maxMember["blip_y"],maxMember["blip_z"]
                nRun,nSubRun,nEvent = maxMember["run"],maxMember["timestamp"],maxMember["event"]

                completeDetDF_nearby_blips = completeDetDF.query('(run == %s) and (timestamp == %s) and(event == %s) and (blip_x > %s and blip_x < %s) and (blip_y > %s and blip_y < %s) and (blip_z>%s and blip_z< %s)'%(nRun, nSubRun, nEvent,blip_x-20.0,blip_x+20.0,blip_y-20.0,blip_y+20.0,blip_z-20.0,blip_z+20.0))
                completeNearbyEs = completeDetDF_nearby_blips["blip_energy"].agg(list)

                nearbyESums = curr_Max_E
                for nearE in completeNearbyEs:
                    if nearE not in usedEnergies:
                        usedEnergies.append(nearE)
                        nearbyESums += nearE
                        multiplicity += 1
                        if nearE in blip_Es_ievent:
                            blip_Es_ievent.remove(nearE)
                multiplicityInRadius.append(multiplicity)
                summedBlipEs.append(nearbyESums)
            elif maxesUsed == numMaxPts:
                break

        counter+=1
        if counter > 1600:
                break
    print "Number of energy clusters in this DF: ", len(summedBlipEs), " out of "," energies total"
    return summedBlipEs,multiplicityInRadius

def returnMinBlipEs(thisDF):  #return minimum blip energies for each event
    minblips = []
    df_grouped = thisDF.groupby(["event","run","timestamp"])
    print "Number of min blip points: ", df_grouped["blip_energy"].min().count() #counting number of unique events
    #Here, we're cycling through all the stuff in one event -- since we're on blip branch, this is a collection of blip points
    for group_name, df_group in df_grouped:
        eventMinE = df_group["blip_energy"].min()
        if eventMinE < -10000.0 or eventMinE > 10000.0:
            continue
        minblips.append(eventMinE)
    return minblips

def returnMinBlipQs(thisDF):  #Return minimum blip charges for an event
    minblips = []
    df_grouped = thisDF.groupby(["event","run","timestamp"])
    print "Number of min blip points: ", df_grouped["blip_cluster_charge"].min().count() #counting number of unique events
    #Here, we're cycling through all the stuff in one event -- since we're on blip branch, this is a collection of blip points
    for group_name, df_group in df_grouped:
        eventMinq = df_group["blip_cluster_charge"].min()
        minblips.append(eventMinq)
    return minblips

def returnMinBlipsAndMultiplicityAllData(thisDF,completeDetDF,summingRadius): #sum blips within summing radius -- in a region, but then add blips from anywhere nearby in the detector volume
    summedBlipEs = [] #for First blip_Q in an event + nearest neighbor + next nearest neighbor
    eventMins = []
    multiplicityInRadius = []
    df_grouped = thisDF.groupby(["event","run","timestamp"])
    #cd_grouped = completeDetDF(["event","run","timestamp"])
    print "Number of events with energies: ", df_grouped["blip_energy"].first().count() #counting number of unique events
    counter = 0
    for group_name, df_group in df_grouped:
        usedEnergies = []
        blip_Es_ievent = df_group["blip_energy"].agg(list)
        for j in range(len(df_group)):
            if len(blip_Es_ievent) > 0:
                curr_Min_E = min(blip_Es_ievent)
                if curr_Min_E < -10000.0 or curr_Min_E > 10000.0:    #Exclude crazy outliers
                    break
                eventMins.append(curr_Min_E)
                multiplicity = 1
                minIndex = blip_Es_ievent.index(curr_Min_E)
                usedEnergies.append(curr_Min_E)
                blip_Es_ievent.remove(curr_Min_E)
                minMember = df_group.iloc[minIndex]
                blip_x,blip_y,blip_z = minMember["blip_x"],minMember["blip_y"],minMember["blip_z"]
                nRun,nSubRun,nEvent = minMember["run"],minMember["timestamp"],minMember["event"]

                completeDetDF_nearby_blips = completeDetDF.query('(run == %s) and (timestamp == %s) and(event == %s) and (blip_x > %s and blip_x < %s) and (blip_y > %s and blip_y < %s) and (blip_z>%s and blip_z< %s)'%(nRun, nSubRun, nEvent,blip_x-20.0,blip_x+20.0,blip_y-20.0,blip_y+20.0,blip_z-20.0,blip_z+20.0))
                completeNearbyEs = completeDetDF_nearby_blips["blip_energy"].agg(list)

                nearbyESums = curr_Min_E
                for nearE in completeNearbyEs:
                    if nearE not in usedEnergies:
                        usedEnergies.append(nearE)
                        nearbyESums += nearE
                        multiplicity += 1
                        if nearE in blip_Es_ievent:
                            blip_Es_ievent.remove(nearE)
                multiplicityInRadius.append(multiplicity)
                summedBlipEs.append(nearbyESums)

        counter+=1
        if counter > 2700:
                break
    print "Number of energy clusters in this DF: ", len(summedBlipEs), " out of ", len(eventEs), " energies total"
    return eventMins, summedBlipEs,multiplicityInRadius


def returnMultiplicityEsBinned(multiplicities,EsForMultiplicities):
    df_mult = pd.DataFrame({'multiplicities':multiplicities, 'Es':EsForMultiplicities})
    return df_mult


def sumMaxClusterEnergiesWithNeighbors(thisDF): #Adds clusters of blips togetehr for a larger total
    summedCloseClusters = [] #for max blip_Q in an event + nearest neighbor
    maxEs = []
    df_grouped = thisDF.groupby(["event","run","timestamp"])
    print "Number of max blip points: ", df_grouped["blip_energy"].max().count() #counting number of unique events
    #Here, we're cycling through all the stuff in one event -- since we're on blip branch, this is a collection of blip points
    for group_name, df_group in df_grouped:
        #print df_group[group_name]
        #print df_group
        #print "\nDf_group EV min: ",df_group["blip_cluster_charge"].min()

        eventMaxE = df_group["blip_energy"].max()
        #print "EventMaxq: ",eventMaxq
        if eventMaxE > -0.01:
            maxLocation = df_group.loc[df_group["blip_energy"].idxmax()]
        else:
            continue

        blip_Es_ievent = df_group["blip_energy"].agg(list)
        distancesfromMax = distancesFromPt(maxLocation,df_group)
        #print "d's: ",distancesfromMax

        if len(distancesfromMax) != len(blip_Es_ievent):
            print "Length mismatch... in distances and q's"
        else:
            maxEs.append(eventMaxE)
            Eindex = distancesfromMax.index(min(distancesfromMax))
            summedCloseClusters.append(eventMaxE+blip_Es_ievent[Eindex])
            #print "Sum q's: ",summedCloseClusters
    return maxEs,summedCloseClusters

def sumNearbyClusterEnergies(thisDF):
    summedBlipEs = [] #for First blip_Q in an event + nearest neighbor + next nearest neighbor
    eventEs = []
    df_grouped = thisDF.groupby(["event","run","timestamp"])
    print "Number of events with energies: ", df_grouped["blip_energy"].first().count() #counting number of unique events
    counter = 0
    for group_name, df_group in df_grouped:
        used_members = [] #stores indices of group member blips already summed up
        for j in range(len(df_group)):
            if j in used_members:#skips blips if already added
                continue
            else:
                used_members.append(j)
                #print df_group.iloc[j]
                member = df_group.iloc[j]
                blipE = member["blip_energy"]
                eventEs.append(blipE)
                nearbyESums = blipE
                if (j+1) < len(df_group):
                    for k in range((j+1),len(df_group)):
                        if k not in used_members:
                            nextMember = df_group.iloc[k]
                            thisDistance = DistTweenTwoPoints(member,nextMember)
                            if thisDistance > 0 and thisDistance < 30:  #If close to first blip, add the energy in, put on list of already added blips
                                nearBlipE = nextMember["blip_energy"]
                                nearbyESums+= nearBlipE
                                eventEs.append(nearBlipE)
                                used_members.append(k)
                summedBlipEs.append(nearbyESums)
        counter+=1
        if counter > 1000:
                break
    print "Number of energy clusters in this DF: ", len(summedBlipEs), " out of ", len(eventEs), " energies total"
    return eventEs, summedBlipEs

def sumMaxAndNeighborsInRadius(thisDF,summingRadius): #sum blip ernergies within some well-defined 3D summing radius (slow)
    summedBlipEs = [] #for First blip_Q in an event + nearest neighbor + next nearest neighbor
    eventEs = []
    multiplicityInRadius = []
    df_grouped = thisDF.groupby(["event","run","timestamp"])
    print "Number of events with energies: ", df_grouped["blip_energy"].first().count() #counting number of unique events
    counter = 0
    for group_name, df_group in df_grouped:
        blip_Es_ievent = df_group["blip_energy"].agg(list)
        print counter
        for j in range(len(df_group)):
            if len(blip_Es_ievent) > 0:
                curr_Max_E = max(blip_Es_ievent)
                #if curr_Max_E < 0.50:    #Only want E > 0.500 MeV
                #    break
                maxIndex = blip_Es_ievent.index(curr_Max_E)
                blip_Es_ievent.remove(curr_Max_E)
                maxMember = df_group.iloc[maxIndex]
                blip_x,blip_y,blip_z = maxMember["blip_x"],maxMember["blip_y"],maxMember["blip_z"]
                eventEs.append(curr_Max_E)
                nearbyESums = curr_Max_E
                multiplicity = 1
                if (j+1) < len(df_group):
                    for k in range((j+1),len(df_group)):
                        nextMember = df_group.iloc[k]
                        nextE = nextMember["blip_energy"]
                        if nextE in blip_Es_ievent:
                            thisDistance = DistTweenMaxNPoint(blip_x,blip_y,blip_z,nextMember)
                            if thisDistance > 0.0 and thisDistance < summingRadius:  #If close to first blip, add the energy in, put on list of already added blips
                                multiplicity += 1
                                blip_Es_ievent.remove(nextE)
                                nearbyESums+= nextE
                                eventEs.append(nextE)
                multiplicityInRadius.append(multiplicity)
                summedBlipEs.append(nearbyESums)

        counter+=1
        if counter > 600:
                break
    print "Number of energy clusters in this DF: ", len(summedBlipEs), " out of ", len(eventEs), " energies total"
    return eventEs, summedBlipEs,multiplicityInRadius

def sumMaxAndNeighborsInRadiusAllData(thisDF,completeDetDF,summingRadius): #sum blips within summing radius -- in a region, but then add blips from anywhere nearby in the detector volume
    summedBlipEs = [] #for First blip_Q in an event + nearest neighbor + next nearest neighbor
    eventEs = []
    multiplicityInRadius = []
    df_grouped = thisDF.groupby(["event","run","timestamp"])
    #cd_grouped = completeDetDF(["event","run","timestamp"])
    print "Number of events with energies: ", df_grouped["blip_energy"].first().count() #counting number of unique events
    counter = 0
    for group_name, df_group in df_grouped:
        usedEnergies = []
        blip_Es_ievent = df_group["blip_energy"].agg(list)
        for j in range(len(df_group)):
            if len(blip_Es_ievent) > 0:
                curr_Max_E = max(blip_Es_ievent)
                eventEs.append(curr_Max_E)
                multiplicity = 1
                if curr_Max_E < 0.50:    #Only want E > 0.500 MeV
                    break
                maxIndex = blip_Es_ievent.index(curr_Max_E)
                usedEnergies.append(curr_Max_E)
                blip_Es_ievent.remove(curr_Max_E)
                maxMember = df_group.iloc[maxIndex]
                blip_x,blip_y,blip_z = maxMember["blip_x"],maxMember["blip_y"],maxMember["blip_z"]
                nRun,nSubRun,nEvent = maxMember["run"],maxMember["timestamp"],maxMember["event"]

                completeDetDF_nearby_blips = completeDetDF.query('(run == %s) and (timestamp == %s) and(event == %s) and (blip_x > %s and blip_x < %s) and (blip_y > %s and blip_y < %s) and (blip_z>%s and blip_z< %s)'%(nRun, nSubRun, nEvent,blip_x-20.0,blip_x+20.0,blip_y-20.0,blip_y+20.0,blip_z-20.0,blip_z+20.0))
                completeNearbyEs = completeDetDF_nearby_blips["blip_energy"].agg(list)

                nearbyESums = curr_Max_E
                for nearE in completeNearbyEs:
                    if nearE not in usedEnergies:
                        usedEnergies.append(nearE)
                        nearbyESums += nearE
                        multiplicity += 1
                        if nearE in blip_Es_ievent:
                            blip_Es_ievent.remove(nearE)
                            eventEs.append(nearE)
                multiplicityInRadius.append(multiplicity)
                summedBlipEs.append(nearbyESums)

        counter+=1
        if counter > 2700:
                break
    print "Number of energy clusters in this DF: ", len(summedBlipEs), " out of ", len(eventEs), " energies total"
    return eventEs, summedBlipEs,multiplicityInRadius

def sumMaxAndNeighborsInRadiusAllData(thisDF,completeDetDF,summingRadius): #sum blips within summing radius -- in a region, but then add blips from anywhere nearby in the detector volume
    summedBlipEs = [] #for First blip_Q in an event + nearest neighbor + next nearest neighbor
    eventEs = []
    multiplicityInRadius = []

    df_grouped = thisDF.groupby(["event","run","timestamp"])
    #cd_grouped = completeDetDF(["event","run","timestamp"])
    print "Number of events with energies: ", df_grouped["blip_energy"].first().count() #counting number of unique events
    counter = 0
    for group_name, df_group in df_grouped:
        usedEnergies = []
        blip_Es_ievent = df_group["blip_energy"].agg(list)
        for j in range(len(df_group)):
            if len(blip_Es_ievent) > 0:
                curr_Max_E = max(blip_Es_ievent)
                eventEs.append(curr_Max_E)
                if curr_Max_E < 0.50:    #Only want E > 0.500 MeV
                    break
                multiplicity = 1
                maxIndex = blip_Es_ievent.index(curr_Max_E)
                usedEnergies.append(curr_Max_E)
                blip_Es_ievent.remove(curr_Max_E)
                maxMember = df_group.iloc[maxIndex]
                blip_x,blip_y,blip_z = maxMember["blip_x"],maxMember["blip_y"],maxMember["blip_z"]
                nRun,nSubRun,nEvent = maxMember["run"],maxMember["timestamp"],maxMember["event"]

                completeDetDF_nearby_blips = completeDetDF.query('(run == %s) and (timestamp == %s) and(event == %s) and (blip_x > %s and blip_x < %s) and (blip_y > %s and blip_y < %s) and (blip_z>%s and blip_z< %s)'%(nRun, nSubRun, nEvent,blip_x-20.0,blip_x+20.0,blip_y-20.0,blip_y+20.0,blip_z-20.0,blip_z+20.0))
                completeNearbyEs = completeDetDF_nearby_blips["blip_energy"].agg(list)

                nearbyESums = curr_Max_E
                for nearE in completeNearbyEs:
                    if nearE not in usedEnergies:
                        usedEnergies.append(nearE)
                        nearbyESums += nearE
                        multiplicity += 1
                        if nearE in blip_Es_ievent:
                            blip_Es_ievent.remove(nearE)
                            eventEs.append(nearE)
                multiplicityInRadius.append(multiplicity)
                summedBlipEs.append(nearbyESums)

        counter+=1
        if counter > 1600:
                break
    print "Number of energy clusters in this DF: ", len(summedBlipEs), " out of ", len(eventEs), " energies total"
    return eventEs, summedBlipEs,multiplicityInRadius

def sumOnlyMultiBlipsEvents(thisDF,summingRadius): #sum blips within summing radius, and only in events with many blips
    summedBlipEs = [] #for First blip_Q in an event + nearest neighbor + next nearest neighbor
    eventEs = []
    multiplicityInRadius = []
    df_grouped = thisDF.groupby(["event","run","timestamp"])
    print "Number of events with energies: ", df_grouped["blip_energy"].first().count() #counting number of unique events
    counter = 0
    for group_name, df_group in df_grouped:
        blip_Es_ievent = df_group["blip_energy"].agg(list)
        for j in range(len(df_group)):
            if len(blip_Es_ievent) > 0:
                curr_Max_E = max(blip_Es_ievent)
                #if curr_Max_E < 0.50:    #Only want E > 0.500 MeV
                #    break
                maxIndex = blip_Es_ievent.index(curr_Max_E)
                blip_Es_ievent.remove(curr_Max_E)
                maxMember = df_group.iloc[maxIndex]
                blip_x,blip_y,blip_z = maxMember["blip_x"],maxMember["blip_y"],maxMember["blip_z"]
                nearbyESums = curr_Max_E
                multiplicity = 1
                if (j+1) < len(df_group):
                    for k in range((j+1),len(df_group)):
                        nextMember = df_group.iloc[k]
                        nextE = nextMember["blip_energy"]
                        if nextE in blip_Es_ievent:
                            thisDistance = DistTweenMaxNPoint(blip_x,blip_y,blip_z,nextMember)
                            if thisDistance > 0.0 and thisDistance < summingRadius:  #If close to first blip, add the energy in, put on list of already added blips
                                multiplicity += 1
                                blip_Es_ievent.remove(nextE)
                                nearbyESums+= nextE
                                eventEs.append(nextE)

                multiplicityInRadius.append(multiplicity)
                if multiplicity > 1:
                    summedBlipEs.append(nearbyESums)

        counter+=1
        if counter > 600:
                break
    print "Number of energy clusters in this DF: ", len(summedBlipEs), " out of ", len(eventEs), " energies total"
    return summedBlipEs,multiplicityInRadius
