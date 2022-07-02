# AZH_studies

- l3KinProducer.py - this version should work with UL sample

## main branch 
- l3KinProducerMAINBranchVersion.py - this was the version that worked earlier on the one file I got from Sahithi - based on the main branch
- l3KinProducerJune21_ValidationNeutrinoSolutionMethod.py - was used to validate neutrino stuff - but was based on the main branch l3Kin
- l3Kin_getGeneratorValues.py - to get generator values - but also it was based on the main branch

## To do 
- need to run over larger sample - get the new generator values, will also have to fix the scripts to work with UL sample
- need scale factors for ttbar sample

# other
- how many neutrino solutions were 'real' as opposed to complex
- ideally, implement ML for neutrino solution

## notes re: changes in branches
- UL now doesn't save the btag for CleanJets. So, need to create a new 4 vector for that
- new version of l3Kin requires an import of itertools, TMath 
