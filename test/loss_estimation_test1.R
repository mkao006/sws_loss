suppressMessages({
    library(faosws)
    library(faoswsUtil)
    library(faoswsProductionImputation)
    library(zoo)
    library(lattice)
    library(data.table)
    library(magrittr)
    library(reshape2)
    library(igraph)
    library(lme4)
})

verbose = FALSE

if(verbose){
    startingTime = Sys.time()
    currentTime = startingTime
}


## Year should be a paramameter selected.
requiredCountries = "840"
selectedYear = as.character(1992:2015)

areaVar = "geographicAreaM49"
yearVar = "timePointYears"
itemVar = "measuredItemCPC"
elementVar = "measuredElement"
valuePrefix = "Value_"
flagObsPrefix = "flagObservationStatus_"
flagMethodPrefix = "flagMethod_"


## Set up testing environments
if(Sys.getenv("USER") == "mk"){
    GetTestEnvironment(
        ## baseUrl = "https://hqlqasws1.hq.un.fao.org:8181/sws",
        baseUrl = "https://hqlprswsas1.hq.un.fao.org:8181/sws",
        token = "09aeec5c-2cf5-46f7-90e7-f8b4bdc6dc93"
        )
    verbose = TRUE
    R_SWS_SHARE_PATH = getwd()
} else {
    R_SWS_SHARE_PATH = "/work/SWS_R_Share/kao"
}



## Function to obtain all CPC item 
getAllItemCPC = function(){
    itemEdgeList =
        adjacent2edge(
            GetCodeTree(domain = "agriculture",
                        dataset = "agriculture",
                        dimension = itemVar)
        )
    itemEdgeGraph = graph.data.frame(itemEdgeList)
    itemDist = shortest.paths(itemEdgeGraph, v = "0", mode = "out")
    fbsItemCodes = colnames(itemDist)[is.finite(itemDist)]
    fbsItemCodes
}

requiredItems = getAllItemCPC()

getProductionData = function(){
    allCountries =
        GetCodeList(domain = "agriculture",
                    dataset = "agriculture",
                    dimension = "geographicAreaM49")[type == "country", code]
    
    productionKey = DatasetKey(
        domain = "agriculture",
        dataset = "agriculture",
        dimensions = list(
            Dimension(name = areaVar,
                      keys = allCountries),
            Dimension(name = elementVar,
                      keys = "5510"),
            Dimension(name = itemVar,
                      keys = requiredItems),
            Dimension(name = yearVar,
                      keys = selectedYear)
        )
    )

    ## Pivot to vectorize yield computation
    productionPivot = c(
        Pivoting(code = areaVar, ascending = TRUE),
        Pivoting(code = itemVar, ascending = TRUE),
        Pivoting(code = yearVar, ascending = FALSE),
        Pivoting(code = elementVar, ascending = TRUE)
    )

    ## Query the data
    productionQuery = GetData(
        key = productionKey,
        flags = TRUE,
        normalized = FALSE,
        pivoting = productionPivot
    )

    ## Convert time to numeric
    productionQuery[, timePointYears := as.numeric(timePointYears)]
    productionQuery

}

getImportData = function(){
    allCountries =
        GetCodeList(domain = "trade",
                    dataset = "total_trade_CPC",
                    dimension = "geographicAreaM49")[type == "country", code]

    importKey = DatasetKey(
        domain = "trade",
        dataset = "total_trade_CPC",
        dimensions = list(
            Dimension(name = areaVar,
                      keys = allCountries),
            Dimension(name = "measuredElementTrade",
                      keys = "5600"),
            Dimension(name = "measuredItemCPC",
                      keys = requiredItems),
            Dimension(name = yearVar,
                      keys = selectedYear)
        )
    )

    ## Pivot to vectorize yield computation
    importPivot = c(
        Pivoting(code = areaVar, ascending = TRUE),
        Pivoting(code = "measuredItemCPC", ascending = TRUE),
        Pivoting(code = yearVar, ascending = FALSE),
        Pivoting(code = "measuredElementTrade", ascending = TRUE)
    )

    ## Query the data
    importQuery = GetData(
        key = importKey,
        flags = TRUE,
        normalized = FALSE,
        pivoting = importPivot
    )

    ## NOTE (Michael): The unit for trade is in kg while for other
    ##                 elements are tonnes, so we divide the trade by
    ##                 1000 to match the unit.
    importQuery[, Value_measuredElementTrade_5600 :=
                   computeRatio(Value_measuredElementTrade_5600, 1000)]
    
    setnames(importQuery,
             old = grep("measuredElementTrade",
                 colnames(importQuery), value = TRUE),
             new = gsub("measuredElementTrade", "measuredElement",
                 grep("measuredElementTrade",
                      colnames(importQuery), value = TRUE)))


    ## Convert time to numeric
    importQuery[, timePointYears := as.numeric(timePointYears)]
    importQuery

}




getOfficialLossData = function(){
    allCountries =
        GetCodeList(domain = "lossWaste",
                    dataset = "loss",
                    dimension = "geographicAreaM49")[type == "country", code]

    ## HACK (Michael): This is a hack, beacause the item hierachy
    ##                 configuration is different in the loss data set
    getLossItemCPC = function(){
        itemEdgeList =
            adjacent2edge(
                GetCodeTree(domain = "lossWaste",
                            dataset = "loss",
                            dimension = "measuredItemSuaFbs")
            )
        itemEdgeGraph = graph.data.frame(itemEdgeList)
        itemDist = shortest.paths(itemEdgeGraph, v = "0", mode = "out")
        fbsItemCodes = colnames(itemDist)[is.finite(itemDist)]
        fbsItemCodes
    }

    lossItems = getLossItemCPC()
   
    
    ## NOTE (Michael): The cpc tree loaded in the loss data set is
    ##                 different to the rest. Thus I can not query
    ##                 item such as 0419.

    lossKey = DatasetKey(
        domain = "lossWaste",
        dataset = "loss",
        dimensions = list(
            Dimension(name = areaVar,
                      keys = allCountries),
            Dimension(name = "measuredElementSuaFbs",
                      keys = "5120"),
            Dimension(name = "measuredItemSuaFbs",
                      keys = lossItems),
            Dimension(name = yearVar,
                      keys = selectedYear)
        )
    )

    ## Pivot to vectorize yield computation
    lossPivot = c(
        Pivoting(code = areaVar, ascending = TRUE),
        Pivoting(code = "measuredItemSuaFbs", ascending = TRUE),
        Pivoting(code = yearVar, ascending = FALSE),
        Pivoting(code = "measuredElementSuaFbs", ascending = TRUE)
    )

    ## Query the data
    lossQuery = GetData(
        key = lossKey,
        flags = TRUE,
        normalized = FALSE,
        pivoting = lossPivot
    )

    setnames(lossQuery,
             old = grep("measuredElementSuaFbs",
                 colnames(lossQuery), value = TRUE),
             new = gsub("measuredElementSuaFbs", "measuredElement",
                 grep("measuredElementSuaFbs",
                      colnames(lossQuery), value = TRUE)))
    setnames(lossQuery,
             old = "measuredItemSuaFbs",
             new = "measuredItemCPC")


    ## Convert time to numeric
    lossQuery[, timePointYears := as.numeric(timePointYears)]
    lossQuery[flagObservationStatus_measuredElement_5120 == "", ]
}


getUnofficialLossData = function(){
    allCountries =
        GetCodeList(domain = "lossWaste",
                    dataset = "loss",
                    dimension = "geographicAreaM49")[type == "country", code]

    ## HACK (Michael): This is a hack, beacause the item hierachy
    ##                 configuration is different in the loss data set
    getLossItemCPC = function(){
        itemEdgeList =
            adjacent2edge(
                GetCodeTree(domain = "lossWaste",
                            dataset = "loss",
                            dimension = "measuredItemSuaFbs")
            )
        itemEdgeGraph = graph.data.frame(itemEdgeList)
        itemDist = shortest.paths(itemEdgeGraph, v = "0", mode = "out")
        fbsItemCodes = colnames(itemDist)[is.finite(itemDist)]
        fbsItemCodes
    }

    lossItems = getLossItemCPC()
   
    
    ## NOTE (Michael): The cpc tree loaded in the loss data set is
    ##                 different to the rest. Thus I can not query
    ##                 item such as 0419.

    lossKey = DatasetKey(
        domain = "lossWaste",
        dataset = "loss",
        dimensions = list(
            Dimension(name = areaVar,
                      keys = allCountries),
            Dimension(name = "measuredElementSuaFbs",
                      keys = "5120"),
            Dimension(name = "measuredItemSuaFbs",
                      keys = lossItems),
            Dimension(name = yearVar,
                      keys = selectedYear)
        )
    )

    ## Pivot to vectorize yield computation
    lossPivot = c(
        Pivoting(code = areaVar, ascending = TRUE),
        Pivoting(code = "measuredItemSuaFbs", ascending = TRUE),
        Pivoting(code = yearVar, ascending = FALSE),
        Pivoting(code = "measuredElementSuaFbs", ascending = TRUE)
    )

    ## Query the data
    lossQuery = GetData(
        key = lossKey,
        flags = TRUE,
        normalized = FALSE,
        pivoting = lossPivot
    )

    setnames(lossQuery,
             old = grep("measuredElementSuaFbs",
                 colnames(lossQuery), value = TRUE),
             new = gsub("measuredElementSuaFbs", "measuredElement",
                 grep("measuredElementSuaFbs",
                      colnames(lossQuery), value = TRUE)))
    setnames(lossQuery,
             old = "measuredItemSuaFbs",
             new = "measuredItemCPC")


    ## Convert time to numeric
    lossQuery[, timePointYears := as.numeric(timePointYears)]
    lossQuery[flagObservationStatus_measuredElement_5120 != "", ]
}

## Need to get national fbs data
getNationalFbs = function(){
    nationalFbs = GetTableData(schemaName = "ess", tableName = "national_fbs")
    setnames(nationalFbs, old = colnames(nationalFbs),
             new = c("geographicAreaM49", "timePointYears",
                 "Value_measuredElement_5510", "Value_measuredElement_5600",
                 "Value_measuredElement_5910", "Value_measuredElement_5712",
                 "Value_measuredElement_5120", "Value_measuredElement_5525",
                 "measuredItemCPC"))
    nationalFbs[, timePointYears := as.numeric(timePointYears)]
    nationalFbs
}

getLossWorldBankData = function(){
    allCountries =
        GetCodeList(domain = "WorldBank",
                    dataset = "wb_ecogrw",
                    dimension = "geographicAreaM49")[type == "country", code]
   
    infrastructureKey =
        DatasetKey(domain = "WorldBank",
                   dataset = "wb_infrastructure",
                   dimensions =
                       list(
                           Dimension(name = "geographicAreaM49",
                                     keys = allCountries),
                           Dimension(name = "wbIndicator",
                                     keys = "IS.ROD.PAVE.ZS"),
                           Dimension(name = "timePointYears",
                                     keys = selectedYear)
                       )
                   )

    climateKey =
        DatasetKey(domain = "WorldBank",
                   dataset = "wb_climate",
                   dimensions =
                       list(
                           Dimension(name = "geographicAreaM49",
                                     keys = allCountries),
                           Dimension(name = "wbIndicator",
                                     keys = c("SWS.FAO.PREC",
                                         "SWS.FAO.TEMP")),
                           Dimension(name = "timePointYears",
                                     keys = selectedYear)
                       )
                   )
    
    gdpKey =
        DatasetKey(domain = "WorldBank",
                   dataset = "wb_ecogrw",
                   dimensions =
                       list(
                           Dimension(name = "geographicAreaM49",
                                     keys = allCountries),
                           Dimension(name = "wbIndicator",
                                     keys = c("NY.GDP.MKTP.PP.KD",
                                         "NY.GDP.PCAP.KD")),
                           Dimension(name = "timePointYears",
                                     keys = selectedYear)
                       )
                   )

    newPivot = c(
        Pivoting(code = "geographicAreaM49", ascending = TRUE),
        Pivoting(code = "wbIndicator", ascending = TRUE),
        Pivoting(code = "timePointYears", ascending = FALSE)
    )

    base =
        data.table(geographicAreaM49 = character(),
                   wbIndicator = character(),
                   timePointYears = character(),
                   Value = numeric())

    merged =
        Reduce(f = function(base, key){
            rbind(base, GetData(key, pivoting = newPivot))
        }, x = list(climateKey, infrastructureKey, gdpKey), init = base)
    
    casted =
        dcast.data.table(merged,
                         geographicAreaM49 + timePointYears ~ wbIndicator,
                         value.var = "Value")
    setnames(casted,
             old = c("IS.ROD.PAVE.ZS", "NY.GDP.MKTP.PP.KD",
                 "NY.GDP.PCAP.KD", "SWS.FAO.PREC", "SWS.FAO.TEMP"),
             new = c("sharePavedRoad", "gdpPPP", "gdpPerCapita",
                 "precipitation", "temperature"))
    casted[, timePointYears := as.numeric(timePointYears)]
    setkeyv(casted, cols = c("geographicAreaM49", "timePointYears"))
    casted
}


## Function to load the loss food group classification
getLossFoodGroup = function(){
    lossFoodGroup = GetTableData(schemaName = "ess", tableName = "loss_food_group")
    setnames(lossFoodGroup, old = colnames(lossFoodGroup),
             new = c("measuredItemFS", "measuredItemNameFS", "foodGroupName",
                 "foodGroup", "foodGeneralGroup", "foodPerishableGroup",
                 "measuredItemCPC"))
    lossFoodGroup[, list(measuredItemCPC, foodGroupName,
                         foodGroup, foodGeneralGroup, foodPerishableGroup)]
    lossFoodGroup
}

## Function to load the loss region classification
getLossRegionClass = function(){
    regionMapping =
        GetTableData(schemaName = "ess", tableName = "loss_region_mapping")
    setnames(regionMapping, old = colnames(regionMapping),
             new = c("geographicAreaM49", "lossRegionClass"))    
    regionMapping
}


imputeSharePavedRoad = function(wbData, pavedRoadVar){
    foo = function(x){
        if(length(na.omit(x)) >= 2){
            tmp = na.locf(na.approx(x, na.rm = FALSE), na.rm = FALSE)
        } else {
            tmp = x
        }
        tmp
    }
    wbData[, `:=`(c(pavedRoadVar),
                  foo(.SD[[pavedRoadVar]])),
         by = "geographicAreaM49"]
}

mergeAllLossData = function(lossData, ...){
    explanatoryData = list(...)
    Reduce(f = function(x, y){
        keys = intersect(colnames(x), colnames(y))
        setkeyv(x, keys)
        setkeyv(y, keys)
        merge(x, y, all.x = TRUE)
    },
           x = explanatoryData, init = lossData
           )
}

production = getProductionData()
import = getImportData()
loss = getOfficialLossData()
## nationalFbs = getNationalFbs()
wb = getLossWorldBankData() %>%
    imputeSharePavedRoad(wbData = ., pavedRoadVar = "sharePavedRoad")
lossFoodGroup = getLossFoodGroup()
lossRegionClass = getLossRegionClass()

finalLoss = copy(loss)

modelData =
    mergeAllLossData(lossData = finalLoss,
                     production, import, wb, lossFoodGroup, lossRegionClass)
setnames(modelData,
         old = c("Value_measuredElement_5510", "Value_measuredElement_5600",
             "Value_measuredElement_5120"),
         new = c("production", "import", "loss"))


## Only take positive loss and primary commodity
modelData = modelData[loss > 0, ]
modelData = modelData[foodGeneralGroup == "primary", ]

## Remove carried forward value, as it reduces correlation
modelData[, variance := var(loss, na.rm = TRUE),
          by = c("geographicAreaM49", "measuredItemCPC")]
modelData[, duplicateValue := duplicated(loss),
          by = c("geographicAreaM49", "measuredItemCPC")]
modelData = modelData[!(variance == 0 & duplicateValue), ]

countryTable =
    GetCodeList(domain = "agriculture",
                dataset = "agriculture",
                dimension = "geographicAreaM49")[type == "country", list(code, description)]
setnames(countryTable, old = c("code", "description"),
         new = c("geographicAreaM49", "geographicAreaM49Name"))


finalModelData =
    merge(modelData, countryTable, by = "geographicAreaM49", all.x = TRUE)
finalModelData[, `:=`(c("duplicateValue", "variance", "temperature",
                        "precipitation", "foodGroup", "foodGeneralGroup",
                        "measuredItemFS","measuredItemNameFS",
                        grep("flag", colnames(finalModelData), value = TRUE)),
                      NULL)]

## Convert keys to factor
finalModelData[, `:=`(c("measuredItemCPC", "foodGroupName", 
                   "foodPerishableGroup", "lossRegionClass", "geographicAreaM49",
                   "geographicAreaM49Name"),
                 lapply(c("measuredItemCPC", "foodGroupName", 
                          "foodPerishableGroup", "lossRegionClass",
                          "geographicAreaM49", "geographicAreaM49Name"),
                        FUN = function(x) as.factor(get(x))
                        )
                 )
          ]

## Remove zero production
finalModelData = finalModelData[production != 0, ]

## Just a test

trainIndex = sample(NROW(finalModelData), NROW(finalModelData) * 0.75)
finalTrainData =
    finalModelData[trainIndex, ]

finalTestData =
    finalModelData[-trainIndex, ]

## Exploratory plots
xyplot(log(loss) ~ log(production)|geographicAreaM49Name,
       data = finalModelData,
       type = c("p", "r"))

xyplot(log(loss) ~ log(production)^2|geographicAreaM49Name,
       data = finalModelData,
       type = c("p", "r"))


xyplot(log(loss) ~ gdpPerCapita|geographicAreaM49Name,
       data = finalModelData,
       type = c("p", "r"))


xyplot(log(loss) ~ log(gdpPerCapita),
       data = finalModelData,
       type = c("p", "r"))


xyplot(log(loss) ~ sharePavedRoad|geographicAreaM49Name,
       data = finalModelData,
       type = c("p", "r"))



## Temporary model (without trade)
## ---------------------------------------------------------------------

## The linear regression model
## Import is excluded here in the model as there are no data
lossLmModel =
    lm(log(loss) ~ log(production + 1) + 
           timePointYears + foodGroupName:geographicAreaM49Name,
    data = finalTrainData)

summary(lossLmModel)


## The linear mixed model

## This is the current model
system.time({
    lossLmeModel =
        lmer(log(loss) ~ timePointYears + log(production) + 
                 (-1 + log(production)|foodPerishableGroup/foodGroupName/measuredItemCPC/geographicAreaM49Name),
             data = finalTrainData)
})


test = finalTrainData[, list(lloss = log(loss), foodGroupName, foodPerishableGroup,
                           lossRegionClass, geographicAreaM49, measuredItemCPC,
                           interact = factor(as.numeric(geographicAreaM49) *
                               as.numeric(measuredItemCPC)))]

summary(lossLmeModel)

## Comparison plot
par(mfrow = c(1, 2))
finalTrainData[, predictedLM := exp(predict(lossLmModel, .SD))]
with(finalTrainData, plot(log(loss), log(predictedLM),
                          xlim = c(0, 15), ylim = c(0, 15)))
abline(a = 0, b = 1, col = "red", lty = 2)

finalTrainData[, predictedLME :=
                   exp(predict(lossLmeModel, .SD, allow.new.levels = TRUE))]
with(finalTrainData, plot(log(loss), log(predictedLME),
                          xlim = c(0, 15), ylim = c(0, 15)))
abline(a = 0, b = 1, col = "red", lty = 2)





## Calculation of R-squared
1 - with(finalTrainData, sum((log(predictedLM) - log(loss))^2))/
    with(finalTrainData, sum((log(loss) - mean(log(loss)))^2))

1 - with(finalTrainData, sum((log(predictedLME) - log(loss))^2))/
    with(finalTrainData, sum((log(loss) - mean(log(loss)))^2))

## Plot the random effects
##
## NOTE (Michael): The order make sense, with animal products with the
##                 lowest loss rate followed by non perishable and no
##                 difference between fruits and semi perishable.
randomEffect = ranef(lossLmeModel, condVar = TRUE)
dotplot(randomEffect)[[4]]


## Plot the fit on the original scale
with(finalTrainData,
     plot(loss, predictedLM, xlim = c(0, exp(15)), ylim = c(0, exp(15))))
abline(a = 0, b = 1, col = "red", lty = 2)
with(finalTrainData,
     plot(loss, predictedLME, xlim = c(0, exp(15)), ylim = c(0, exp(15))))
abline(a = 0, b = 1, col = "red", lty = 2)

## Predict out of sample
lmNewLevelPredict = function(model, row){
    tmp = predict(model, newdata = row)
    if(inherits(tmp, "try-error")){
        tmp = NA
    }
    tmp
}

for(i in 1:NROW(finalTestData)){
    finalTestData[i, predictedLM :=
                      exp(lmNewLevelPredict(lossLmModel, row = .SD))]
}
    
finalTestData[production!= 0, predictedLM := exp(predict(lossLmModel, .SD))]
with(finalTestData, plot(log(loss), log(predictedLM), xlim = c(0, 15),
                         ylim = c(0, 15)))
abline(a = 0, b = 1, col = "red", lty = 2)

finalTestData[production != 0, predictedLME :=
                  exp(predict(lossLmeModel, .SD, , allow.new.levels = TRUE))]
with(finalTestData, plot(log(loss), log(predictedLME), xlim = c(0, 15),
                         ylim = c(0, 15)))
abline(a = 0, b = 1, col = "red", lty = 2)


## normal scale comparison
with(finalTestData, plot(loss, predictedLM, xlim = c(0, exp(15)),
                         ylim = c(0, exp(15))))
abline(a = 0, b = 1, col = "red", lty = 2)
with(finalTestData, plot(loss, predictedLME, xlim = c(0, exp(15)),
                         ylim = c(0, exp(15))))
abline(a = 0, b = 1, col = "red", lty = 2)




## Final Model (when trade is completed)
## ---------------------------------------------------------------------

## A model just for 2010, this is the final model when we have trade data

finalModelData2010 = finalModelData[timePointYears == 2010, ]
## finalModelData2010 =
##     finalModelData2010[!geographicAreaM49 %in% c("196", "300", "470", "620"), ]
lmData = finalModelData2010[production != 0 & foodGroupName != "eggs", ]
## Import is excluded in the model as there are no data, time is very
## insignificant.
lossLmModel =
    lm(log(loss) ~ log(production + 1) + import + timePointYears +
           foodGroupName:geographicAreaM49,
    data = lmData)
summary(lossLmModel)

lmData[, predicted := exp(predict(lossLmModel, lmData))]
with(lmData, plot(log(loss), log(predicted), xlim = c(0, 15), ylim = c(0, 15)))
abline(a = 0, b = 1, col = "red", lty = 2)


lmeData =
    finalModelData2010[production != 0, list(geographicAreaM49,
                           geographicAreaM49Name, import, lossRegionClass,
                       timePointYears, measuredItemCPC, production, gdpPerCapita,
                       sharePavedRoad, foodPerishableGroup, foodGroupName, loss)]

## Can not be estimated because we don't have enough data.
lossLmeModel =
    lmer(log(loss) ~ timePointYears + log(import) + log(production) + 
         (-1 + log(production)|foodPerishableGroup/foodGroupName/measuredItemCPC/geographicAreaM49Name),
         data = lmeData)

summary(lossLmeModel)

lmeData[, predicted := exp(predict(lossLmeModel, .SD, allow.new.levels = TRUE))]
with(lmeData, plot(log(loss), log(predicted), xlim = c(0, 15), ylim = c(0, 15)))
abline(a = 0, b = 1, col = "red", lty = 2)

randomEffect = ranef(lossLmeModel, condVar = TRUE)







## Train data
lmMSE = lmData[!is.na(predicted), sqrt(sum((loss - predicted)^2))/.N,
    by = "measuredItemCPC"]
setnames(lmMSE, "V1", "lm")
lmeMSE = lmeData[!is.na(predicted), sqrt(sum((loss - predicted)^2))/.N,
    by = "measuredItemCPC"]
setnames(lmeMSE, "V1", "lme")

compareMSE = merge(lmMSE, lmeMSE, by = "measuredItemCPC")
compareMSE[measuredItemCPC != "0111", list(lm = sum(lm), lme = sum(lme))]
with(compareMSE, plot(log(lm), log(lme)))
abline(0, 1, col = "red", lty = 2)



## For test data
lmMSE = finalTestData[!is.na(predictedLM), sqrt(sum((loss - predictedLM)^2))/.N,
    by = "measuredItemCPC"]
setnames(lmMSE, "V1", "lm")
lmeMSE = finalTestData[!is.na(predictedLME), sqrt(sum((loss - predictedLME)^2))/.N,
    by = "measuredItemCPC"]
setnames(lmeMSE, "V1", "lme")

compareMSE = merge(lmMSE, lmeMSE, by = "measuredItemCPC")
compareMSE[measuredItemCPC != "0111", list(lm = sum(lm), lme = sum(lme))]
with(compareMSE, plot(log(lm), log(lme)))
abline(0, 1, col = "red", lty = 2)

## The result shows that LME is probably more robust as expected.
library(car)
scatterplotMatrix(~loss + predictedLM + predictedLME,
                  data = finalTestData[predictedLM != predictedLME &
                                       geographicAreaM49 != "724", ])

scatterplotMatrix(~loss + predictedLM + predictedLME,
                  data = finalTestData[geographicAreaM49 != "724", ])


check =
    finalTestData[predictedLM != predictedLME,
                  list(loss, predictedLM, predictedLME)]
check[, lmpctdiff := (predictedLM - loss)/loss]
check[, lmepctdiff := (predictedLME - loss)/loss]
print(check[order(lmpctdiff), ], nrow = 1603)


dotplot(~log(lm) + log(lme)|measuredItemCPC, data = compareMSE, auto.key = TRUE)


## Perishable group is an aggregated group of group names, so don't
## use both together.


## Try to add in share of gdp of agriculture

getAllLossData =function(){
    allCountries =
        GetCodeList(domain = "lossWaste",
                    dataset = "loss",
                    dimension = "geographicAreaM49")[type == "country", code]

    ## HACK (Michael): This is a hack, beacause the item hierachy
    ##                 configuration is different in the loss data set
    getLossItemCPC = function(){
        itemEdgeList =
            adjacent2edge(
                GetCodeTree(domain = "lossWaste",
                            dataset = "loss",
                            dimension = "measuredItemSuaFbs")
            )
        itemEdgeGraph = graph.data.frame(itemEdgeList)
        itemDist = shortest.paths(itemEdgeGraph, v = "0", mode = "out")
        fbsItemCodes = colnames(itemDist)[is.finite(itemDist)]
        fbsItemCodes
    }

    lossItems = getLossItemCPC()
   
    
    ## NOTE (Michael): The cpc tree loaded in the loss data set is
    ##                 different to the rest. Thus I can not query
    ##                 item such as 0419.

    lossKey = DatasetKey(
        domain = "lossWaste",
        dataset = "loss",
        dimensions = list(
            Dimension(name = areaVar,
                      keys = allCountries),
            Dimension(name = "measuredElementSuaFbs",
                      keys = "5120"),
            Dimension(name = "measuredItemSuaFbs",
                      keys = lossItems),
            Dimension(name = yearVar,
                      keys = selectedYear)
        )
    )

    ## Pivot to vectorize yield computation
    lossPivot = c(
        Pivoting(code = areaVar, ascending = TRUE),
        Pivoting(code = "measuredItemSuaFbs", ascending = TRUE),
        Pivoting(code = yearVar, ascending = FALSE),
        Pivoting(code = "measuredElementSuaFbs", ascending = TRUE)
    )

    ## Query the data
    lossQuery = GetData(
        key = lossKey,
        flags = TRUE,
        normalized = FALSE,
        pivoting = lossPivot
    )

    setnames(lossQuery,
             old = grep("measuredElementSuaFbs",
                 colnames(lossQuery), value = TRUE),
             new = gsub("measuredElementSuaFbs", "measuredElement",
                 grep("measuredElementSuaFbs",
                      colnames(lossQuery), value = TRUE)))
    setnames(lossQuery,
             old = "measuredItemSuaFbs",
             new = "measuredItemCPC")


    ## Convert time to numeric
    lossQuery[, timePointYears := as.numeric(timePointYears)]
    lossQuery
}

production = getProductionData()
import = getImportData()
loss = getAllLossData()
checkRates = mergeAllLossData(lossData = loss, production, import)

checkRates2010 = checkRates[timePointYears == 2010, ]
checkRates2010[is.na(Value_measuredElement_5510),
               Value_measuredElement_5510 := 0]
checkRates2010[is.na(Value_measuredElement_5600),
               Value_measuredElement_5600 := 0]

officialRates2010 =
    checkRates2010[flagObservationStatus_measuredElement_5120 == "",
                   list(official = mean(computeRatio(Value_measuredElement_5120,
                    (Value_measuredElement_5510 + Value_measuredElement_5600)) * 100,
               na.rm = TRUE)), by = c("measuredItemCPC")]


unofficialRates2010 =
    checkRates2010[flagObservationStatus_measuredElement_5120 != "",
               list(unofficial = mean(computeRatio(Value_measuredElement_5120, 
                    (Value_measuredElement_5510 + Value_measuredElement_5600)) * 100,
               na.rm = TRUE)), by = c("measuredItemCPC")]



ratesComparison =
    merge(officialRates2010, unofficialRates2010, by = "measuredItemCPC")

with(ratesComparison[unofficial <= 100, ],
     plot(official, unofficial, xlim = c(0, 0.5), ylim = c(0, 0.5)))

ratesComparison[measuredItemCPC == "0111", ]
checkRates2010[measuredItemCPC == "0111" &
                 flagObservationStatus_measuredElement_5120 == "", ]
checkRates2010[measuredItemCPC == "0111" &
                 flagObservationStatus_measuredElement_5120 == "", ]


checkRates2010[measuredItemCPC == "01341" &
                 flagObservationStatus_measuredElement_5120 == "", ]
