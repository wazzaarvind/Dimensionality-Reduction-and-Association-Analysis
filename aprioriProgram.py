import numpy as np
from collections import Counter
from itertools import *
import itertools
import time 

itemCount=dict()
copyGeneList=[]
arr=[] #List of list of sets Full Database. Basically a 2D matrix where every gene is a set. 
listSetDatabase=[]
freqSetList = []
ruleBody = []
ruleHead = []
ruleIndex=[]

def mainmethod():
    totalCount=0
    f=open('associationruletestdata.txt','r')

    i=0
    count=0

    #Converting "Up Down Down Up" to "G1_Up G2_Down G3_Down G4_Up"
    for line in f:
        count=1
        setDatabase=set()
        s=[x for x in line.strip().split("\t")]
        arr.append([])
        for word in s:
            if count==len(s):
                arr[i].append(set([word]))
                setDatabase.add(word)
            else:
                updatedWord="G"+str(count)+"_"+word
                setDatabase.add(updatedWord)
                arr[i].append(set([updatedWord]))
                count+=1
        listSetDatabase.append(setDatabase)
        i+=1

    geneList=[]
    support=[0.6]
    print("Support is set to be "+str(int(support[0]*100))+"%")
    minSupport=[]

    #Calculate min support for all input
    for unitSupport in support:
        minSupport.append(unitSupport*len(arr))


    # For Single Element sets alone
    for pos in range(0,len(arr)):
        for oneGeneSet in arr[pos]:
            oneGeneString = ','.join(str(setElement) for setElement in oneGeneSet)
            if oneGeneString in itemCount:
                itemCount[oneGeneString]+=1
            else:
                itemCount[oneGeneString]=1    
    
    listOfGenes=[]
    tempCount = 0

    everyTimeDict = dict()

    #Looping for 1 to 2 alone
    for oneGeneStringItem in itemCount:
        if(itemCount[oneGeneStringItem]>=minSupport[0]):
            geneList.append(oneGeneStringItem)
    n=1
    totalCount+=len(geneList) 
    print("number of length-1 frequent itemsets:",len(geneList))
    freqSetList.append(geneList)


    n+=1

    for i in range(0,len(geneList)):
        for j in range(i+1,len(geneList)):
            everyTimeDict[geneList[i]+","+geneList[j]]=0
    updateDictionary(everyTimeDict)

    #Loop block for generating 3 to n genes in a set
    for supportValue in minSupport:
        while(1):
            geneList=[]
            for oneGeneDictString in itemCount:
                #Find the length of the gene fragments
                listOfGenes=[]
                listOfGenes=oneGeneDictString.split(",") #Finds out the number of genes in oneGeneSet
                if len(listOfGenes)==n:    #If the length of listOfGenes is the same as n of some iteration then we proceed
                    if itemCount[oneGeneDictString]>=supportValue: #If the dictionary contains a value greater than the supportValue then we need to:
                        geneList.append(oneGeneDictString) #add the key or oneGeneSet to the geneList, Here oneGeneSet is a String

            #Exit Condition after checking the length of the geneList
            totalCount+=len(geneList)
            print("number of length-"+str(n)+" frequent itemsets: "+str(len(geneList)))
            if(len(geneList)==0):
                print("number of all length frequent itemsets:",totalCount)
                return
            freqSetList.append(geneList)
            #Incase of no exit
            everyTimeDict = dict()
            #Generate combinations from the list of genesets in the geneList
            for i in range(0,len(geneList)):
                set1 = set(geneList[i].split(","))
                for j in range(i+1,len(geneList)):
                    set2 = set(geneList[j].split(","))
                    intersectedSet = set1&set2
                    if(len(intersectedSet)>=n-1):
                        unionSet=set1|set2
                        innergeneList = list(unionSet)
                        innergeneList.sort()
                        currentGeneString = ','.join(str(setElement) for setElement in innergeneList)
                        everyTimeDict[currentGeneString] = 0

            updateDictionary(everyTimeDict)
            n+=1
                

def updateDictionary(listSetGenes): #List of set of genes #listSetGenes will now be a dictionary
    #When called the first time it is a list of list of sets whereas the second time it becomes a list of sets. 
    #If the iterationList build code is made a function it can again be called seperately 
    for oneSet in listSetGenes: #oneGene is a set and oneRow is a list of set 
        setToCompare = set(oneSet.split(","))
        for pos in range(0,len(listSetDatabase)):
           if setToCompare < listSetDatabase[pos]:
                listToIterate = list(setToCompare)
                listToIterate.sort()
                geneString = ','.join(str(setElement) for setElement in listToIterate)
                if geneString in itemCount:
                    itemCount[geneString]+=1
                else:
                    itemCount[geneString]=1
    return


def templateOne(query):
    query = query.split('(', 1)[1]
    partsOfQuery = []
    partsOfQuery = query.split(",")
    sizeOfParts = len(partsOfQuery)
    listOfQueryGenes = []
    setOfQueryGenes = set()
    partsOfQuery[2] = partsOfQuery[2].split('[')[1]
    partsOfQuery[len(partsOfQuery)-1] = partsOfQuery[len(partsOfQuery)-1].split(']')[0]
    listOfQueryGenes = []
    for i in range(2,len(partsOfQuery)):
        partsOfQuery[i] = partsOfQuery[i].replace("'","")
        listOfQueryGenes.append([partsOfQuery[i]])

    count = 0
    innerCount=0
    for pos in range(0,len(ruleBody)):
        innerCount=0
        for queryGene in listOfQueryGenes:
            bodyList=[]
            headList=[]
            bodyList=ruleBody[pos].split(",")
            headList=ruleHead[pos].split(",")

            #    BLOCK ANY
            if(partsOfQuery[1]==' "ANY"'):
                if partsOfQuery[0]=='"RULE"':
                    if set(queryGene) <= set(bodyList) or set(queryGene) <= set(headList):
                        count = count+1
                        ruleIndex.append(pos)
                        break

                elif partsOfQuery[0]=='"BODY"':
                    if set(queryGene) <= set(bodyList):
                        count = count+1
                        ruleIndex.append(pos)
                        break

                elif partsOfQuery[0]=='"HEAD"':
                    if set(queryGene) <= set(headList):
                        count = count+1
                        ruleIndex.append(pos)
                        break

            

            #    BLOCK NONE
            elif(partsOfQuery[1]==' "NONE"'):
                if partsOfQuery[0]=='"RULE"':
                    if not(set(queryGene) <= set(bodyList)) and not(set(queryGene) <= set(headList)):
                        innerCount+=1
                        if(innerCount==len(listOfQueryGenes)):
                            count+=1
                            ruleIndex.append(pos)

                elif partsOfQuery[0]=='"BODY"':
                    if not(set(queryGene) <= set(bodyList)):
                        innerCount+=1
                        if(innerCount==len(listOfQueryGenes)):
                            count+=1
                            ruleIndex.append(pos)
                
                elif partsOfQuery[0]=='"HEAD"':
                    if not(set(queryGene) <= set(headList)):
                        innerCount+=1
                        if(innerCount==len(listOfQueryGenes)):
                            count+=1
                            ruleIndex.append(pos)



            #    BLOCK 1
            elif(partsOfQuery[1]==' 1'):
                if partsOfQuery[0]=='"RULE"':
                    if set(queryGene) <= set(bodyList) or set(queryGene) <= set(headList):
                        innerCount+=1
                        if innerCount>1:
                            break
                    
                    if innerCount==1 and queryGene==listOfQueryGenes[len(listOfQueryGenes)-1]:
                        count+=1
                        ruleIndex.append(pos)

                elif partsOfQuery[0]=='"BODY"':
                    if set(queryGene) <= set(bodyList):
                        innerCount+=1
                        if innerCount>1:
                            break
                    
                    if innerCount==1 and queryGene==listOfQueryGenes[len(listOfQueryGenes)-1]:
                        count+=1
                        ruleIndex.append(pos)

                elif partsOfQuery[0]=='"HEAD"':
                    if set(queryGene) <= set(headList):
                        innerCount+=1
                        if innerCount>1:
                            break
                    
                    if innerCount==1 and queryGene==listOfQueryGenes[len(listOfQueryGenes)-1]:
                        count+=1
                        ruleIndex.append(pos)

    return count       



def templateTwo(query):
    query = query.split('(', 1)[1]
    partsOfQuery = []
    partsOfQuery = query.split(",")
    sizeOfParts = len(partsOfQuery)
    partsOfQuery[1] = partsOfQuery[1].split(')')[0]
    count=0
    for pos in range(0,len(ruleBody)):
        bodyList=[]
        headList=[]
        bodyList=ruleBody[pos].split(",")
        headList=ruleHead[pos].split(",")
        if partsOfQuery[0]=='"RULE"':
            lengthRule=len(bodyList)+len(headList)
            if lengthRule>=int(partsOfQuery[1]):
                count+=1
                ruleIndex.append(pos)
        
        if partsOfQuery[0]=='"HEAD"':
            lengthRule=len(headList)
            if lengthRule>=int(partsOfQuery[1]):
                count+=1
                ruleIndex.append(pos)

        if partsOfQuery[0]=='"BODY"':
            lengthRule=len(bodyList)
            if lengthRule>=int(partsOfQuery[1]):
                count+=1
                ruleIndex.append(pos)

    return count


def ruleGeneration():
    
    for i in range(len(freqSetList),0,-1): #for each level, where each level has length n strings, n=1,2....

        for eachString in freqSetList[i-1]: #considering gene combinations at each level one by one. Considering each string at n level, like 'abc'.
            tempStringList = []
            tempStringList.extend(eachString.split(',')) #making a list out of the genes we consider, split because A,B,C is eachString, which we ABC.DEF etc
            tempStringSet = set(tempStringList) # converting above list into a set
            confidenceNumerator = itemCount[eachString] #all items in set 
            j = i
            globalListOfStrings = []
            while(j>1):
                    tempList = list(itertools.combinations(tempStringList,j-1))
                    globalListOfStrings.append(tempList)
                    j = j-1

            for listOfStrings in globalListOfStrings:
                for eachStringInList in listOfStrings: #considering each n letter string, like 'ABCD'
                    leftSideStringSet = set()
                    for setElement in eachStringInList: #for A or B or C or D
                        leftSideStringSet.add(setElement) #basically converting 'ABCD' to {A,B,C,D}

                    if(i>1):
                        singleString = ','.join(str(setElement) for setElement in eachStringInList)
                    else:
                        singleString = eachStringInList

                    confidenceDenominator = itemCount[singleString]
                    confidence = confidenceNumerator/confidenceDenominator
                    if(confidence>=0.7):
                        ruleBody.append(singleString)
                        tempRuleHead = tempStringSet - set(eachStringInList)
                        tempRuleHeadList = list(tempRuleHead)
                        ruleHead.append(','.join(str(setElement) for setElement in tempRuleHeadList))

    if len(ruleBody)>0:
        print("---------------------------------- Association Rules ----------------------------------")
    for i in range(0,len(ruleBody)):
        print(ruleBody[i]+"--------->"+ruleHead[i])
    print("The total number of rules generated are", len(ruleBody))

    query = input("\nEnter your query :")
    query = query.split('template', 1)[1]
    if query[0]=='1':
        count=templateOne(query)
        for pos in ruleIndex:
            print(ruleBody[pos]," ----> ",ruleHead[pos])    
        print("The number of rules that match the query is",count)
    elif query[0]=='2':
        count=templateTwo(query)
        for pos in ruleIndex:
            print(ruleBody[pos]," ----> ",ruleHead[pos])
        print("The number of rules that match the query is",count)
    elif query[0]=='3':
        count3=0
        query = query.split('(', 1)[1]
        partsOfQuery = []
        partsOfQuery = query.split(",")
        sizeOfParts = len(partsOfQuery)

        
        if partsOfQuery[0][1]=='1':

            query="("+partsOfQuery[1].replace(" ","")+","+partsOfQuery[2]
            for i in range(3,len(partsOfQuery)-2):
                query+=','+partsOfQuery[i]
            query = query.split(']', 1)[0]
            query+=']'
            query+=')'
            print(query)
            count1=templateOne(query)
            if partsOfQuery[0][len(partsOfQuery[0])-2]=='1':
                query="("+partsOfQuery[4].replace(" ","")+","+partsOfQuery[5]
                for i in range(6,len(partsOfQuery)):
                    query+=','+partsOfQuery[i]
                query = query.split('"', 1)[1]
                query = ''.join(('("',query))

                print(query)
                count2=templateOne(query)
            elif partsOfQuery[0][len(partsOfQuery[0])-2]=='2':    
                query="("+partsOfQuery[len(partsOfQuery)-2].replace(" ","")+","+partsOfQuery[len(partsOfQuery)-1]
                print(query)
                count2=templateTwo(query)

        
        elif  partsOfQuery[0][1]=='2':
            query="("+partsOfQuery[1].replace(" ","")+","+partsOfQuery[2]+")"
            print(query)
            count1=templateTwo(query)
            print(count1)
            if partsOfQuery[0][len(partsOfQuery[0])-2]=='1':
                query="("+partsOfQuery[3].replace(" ","")+","+partsOfQuery[4]
                for i in range(5,len(partsOfQuery)):
                    query+=','+partsOfQuery[i]
                print(query)

                count2=templateOne(query)
            elif partsOfQuery[0][len(partsOfQuery[0])-2]=='2':    
                query="("+partsOfQuery[3].replace(" ","")+","+partsOfQuery[4]
                count2=templateTwo(query)
                print(count2)

        if partsOfQuery[0][2]=='a':
            uniqueIndex=dict()
            for pos in ruleIndex:
                if pos in uniqueIndex:
                    uniqueIndex[pos]+=1
                else:
                    uniqueIndex[pos]=1
            for pos in uniqueIndex:
                if uniqueIndex[pos]==2:
                    print(ruleBody[pos]," ----> ",ruleHead[pos])
                    count3+=1
        elif partsOfQuery[0][2]=='o':
            uniqueIndex=dict()
            for pos in ruleIndex:
                if pos in uniqueIndex:
                    uniqueIndex[pos]+=1
                else:
                    uniqueIndex[pos]=1
            for pos in uniqueIndex:
                if uniqueIndex[pos]>=1:
                    print(ruleBody[pos]," ----> ",ruleHead[pos])
                    count3+=1

        print("The number of rules that match the query is",count3)
        


startTime=time.time()
mainmethod()
ruleGeneration()


