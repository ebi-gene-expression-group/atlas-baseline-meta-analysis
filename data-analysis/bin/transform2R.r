#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(XML))

option_list = list(
  make_option(
    c("-c", "--countstsv"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to a tabular file with counts, presumably library size normalised."
  ),
  make_option(
    c("--sdrf"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to SDRF file."
  ),
  make_option(
    c("--configuration"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to Atlas experiment configuration file."
  ),
  make_option(
    c("--output"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to output RDS file."
  ),
  make_option(
    c("--batch"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Batch variable in the metadata"
  ),
  make_option(
    c("--methods"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to analysis methods file."
  ),
  make_option(
    c("--remove_single_batch_genes"),
    action = "store_true",
    default = FALSE,
    help = 'Remove all zero genes and single batch genes.'
  )
)

op<-OptionParser(option_list=option_list)
opt <- parse_args(op)

suppressPackageStartupMessages(require(ExpressionAtlasInternal))
suppressPackageStartupMessages(require(SummarizedExperiment))

remove_genes_with_expression_on_single_batch<-function(experiment) {
  ma<-NULL
  batches<-unique(experiment$batch)
  for ( batch in batches ) {
    rowSums(experiment[,colData(experiment)$batch == batch]@assays$data[[1]])->batchGeneSum
    if(is.null(ma)) {
      ma<-as.matrix(batchGeneSum)
    } else {
      ma<-cbind(ma, batchGeneSum)
    }
  }
  # remove rows that only have zeros, due to
  # https://github.com/brentp/combat.py/issues/11
  # producing all NaNs
  #rowSums(experimentSummary@assays$data[[1]]) == 0 -> zeroRows
  #print(paste0("Removing ", length(which(zeroRows)), " genes with zero expression on all samples."))
  #experimentSummary<-experimentSummary[!zeroRows, ]

  # ma holds on each column the sum of within samples of a batch
  # ma-row > 0 will give a boolean vector with trues when >0
  # sum(ma-row > 0, na.rm = TRUE) will count how many trues we have.
  # We want all the genes that are present in more than a single batch.
  genesMoreOneBatch<-apply(ma, 1, function(x) sum(x > 0, na.rm = TRUE) > 1 )
  print(paste0("Removed ",sum(!genesMoreOneBatch)," genes that were only present in a single batch."))
  return(experiment[genesMoreOneBatch, ])
}


# parseSDRF

#     - Take an SDRF filename and the experiment type from the XML config, and

#     return a subset of the SDRF containing only the assay names (or ENA runs),

#     the Characteristics, and FactorValue columns. 

#     - It tries to include unit columns as well.

#     - It renames the columns to remove e.g. Characteristics[].

#     - It removes duplicated columns, e.g. if genotype is a Characteristic and a

#     Factor.

#     - It returns the new "SDRF" as a data frame.

parseSDRF <- function( filename, atlasExperimentType ) {



    # Read in the SDRF file. Set header=FALSE because we don't want the column

    # headings to be made "R-safe" -- this confuses things when we're trying to

    # find the Charactersitic and Factor names. Set stringsAsFactors=FALSE so

    # that we can use grep.

    completeSDRF <- read.delim( filename, header = FALSE, stringsAsFactors = FALSE )

    

    cat( "Getting characteristic and factor column indices\n" )



    # Get the Characteristics column indices, and any unit columns next to them.

    charColIndices <- grep( "^Characteristics", ignore.case = FALSE, completeSDRF[ 1, ] )

    charColIndices <- .addUnitCols( charColIndices, completeSDRF )

    

    # Get the Factor column indices, an any unit columns next to them.

    factorColIndices <- grep( "^Factor\\s?Value", ignore.case = FALSE, completeSDRF[ 1, ] )

    factorColIndices <- .addUnitCols( factorColIndices, completeSDRF )



    cat( "Checking for technical replicates\n" )



    # Get the index of the Comment[technical replicate group] column, if there

    # is one.

    techRepGroupColIndex <- grep( "Comment\\s?\\[\\s?technical[ _]replicate[ _]group\\s?\\]", ignore.case = TRUE, completeSDRF[ 1, ] )

    

    cat( "Locating assay names...\n" )



    # Get the column index for assay names. For microarray data, this is "Assay

    # Name" or "Hybridization Name". For RNA-seq data, this is "Comment[ENA_RUN]"

    if( grepl( "rnaseq", atlasExperimentType ) ) {

        

        assayNameColIndex <- grep( "Comment\\s?\\[\\s?ENA_RUN\\s?\\]", ignore.case = TRUE, completeSDRF[ 1, ] )



        if( length( assayNameColIndex ) != 1 ) {

            assayNameColIndex <- grep( "Comment\\s?\\[\\s?RUN_NAME\\s?\\]", ignore.case = TRUE, completeSDRF[ 1, ] )

        }



        if( length( assayNameColIndex ) != 1 ) {

            assayNameColIndex <- grep( "Scan\\s?Name", ignore.case = TRUE, completeSDRF[ 1, ] )

        }

        

        if( length( assayNameColIndex ) != 1 ) {

            stop( "Did not find Comment[ ENA_RUN ] or Comment[ RUN_NAME ] or Scan Name column in SDRF." )

        }

    }

    else {

        

        assayNameColIndex <- grep( "Assay\\s?Name", ignore.case = TRUE, completeSDRF[ 1, ] )



        if( length( assayNameColIndex ) != 1 ) {



            assayNameColIndex <- grep( "Hybridi[sz]ation\\s?Name", ignore.case = TRUE, completeSDRF[ 1, ] )



            # Check that we got something.

            if( length( assayNameColIndex ) != 1 ) {

                stop( "Did not find an Assay Name column or a Hybridization Name column, cannot continue." )

            }

        }



        # For two-colour array data, also want to get the label column.

        if( grepl( "2colour", atlasExperimentType ) ) {

          labelColIndex <- which( completeSDRF[ 1, ] == "Label" )

        }

    }

    

    cat( "Found assay names. Subsetting SDRF...\n" )



    # Now we should have everything we need to get the right columns and make a

    # more friendly SDRF.

    if( grepl( "2colour", atlasExperimentType ) ) {



        subsetSDRF <- completeSDRF[ , c( assayNameColIndex, labelColIndex, charColIndices, factorColIndices ) ]

    }

    else {

        

        subsetSDRF <- completeSDRF[ , c( assayNameColIndex, charColIndices, factorColIndices ) ]

    }



    cat( "Finished subsetting SDRF.\n" )



    # If we got a technical replicate group column, add this at the end of the

    # subsetSDRF.

    if( length( techRepGroupColIndex ) > 0 ) {

      subsetSDRF <- cbind( subsetSDRF, completeSDRF[ , techRepGroupColIndex ] )

    }

                    studyColIndex <- grep("Comment\\s?\\[\\s?Study\\s?\\]", ignore.case = TRUE, completeSDRF[ 1, ] )

# If we got a study group column, add this at the end of the

    # subsetSDRF.

    if( length( studyColIndex ) > 0 ) {

      subsetSDRF <- cbind( subsetSDRF, completeSDRF[ , studyColIndex ] )

    }



    # Next, merge the contents of unit columns with the column before.

    subsetSDRF <- .mergeUnits( subsetSDRF )

    

    cat( "Fixing column headings...\n" )



    # Next thing is to name the columns so they have nice names.

    newColNames <- gsub( "Characteristics\\s?\\[", "", subsetSDRF[1,] )

    newColNames <- gsub( "Factor\\s?Value\\s?\\[", "", newColNames )

    newColNames <- gsub( "\\s?\\]", "", newColNames )

    newColNames[ 1 ] <- "AssayName"



    # Replace the last column name with one for the technical replicate group

    # column, if there is one.

    if( length( techRepGroupColIndex ) > 0 ) {

      newColNames[ length( newColNames ) -1 ] <- "technical_replicate_group"

    }

    
if( length( studyColIndex ) > 0 ) {

     newColNames[ length( newColNames ) ] <- "Study"

    }



    # Replace spaces with underscores.

    newColNames <- gsub( " ", "_", newColNames )

    

    cat( "Finished fixing column headings.\n" )



    cat( "Removing duplicated columns...\n" )



    # Now we've got the new names for the columns, check if any are the same

    # (use "tolower" function to convert all to lower case).

    duplicateColIndices <- which( duplicated( tolower( newColNames ) ) )



    # Remove the duplicated columns from the new column names and the subset

    # SDRF.

    if( length( duplicateColIndices ) > 0 ) {

        subsetSDRF <- subsetSDRF[ , -duplicateColIndices ]

        newColNames <- newColNames[ -duplicateColIndices ]

    }



    cat( "Finished removing duplicated columns.\n" )



    # Remove the first row of the SDRF (this is the old column headings)

    subsetSDRF <- subsetSDRF[ -1, ]

    

    cat( "Applying new column headings...\n" )



    # Add the new column names as the column headings.

    colnames( subsetSDRF ) <- newColNames



    cat( "Finished applying new column headings.\n" )

    

    cat( "Removing duplicated rows...\n" )



    # Remove duplicated rows, which occur e.g. if an assay has more than one file.

    duplicateRowIndices <- which( duplicated( subsetSDRF ) )

    if( length( duplicateRowIndices ) > 0 ) {

        subsetSDRF <- subsetSDRF[ -duplicateRowIndices, ]

    }



    cat( "Finished removing duplicated rows.\n" )



    # Make assay names "R-safe".

    subsetSDRF$AssayName <- make.names( subsetSDRF$AssayName )



    # Return the subset SDRF.

    return( subsetSDRF )

}

# .addUnitCols

#     - Given a vector of column indices and a data frame with the complete SDRF,

#     return a vector of column indices containing the original ones plus any

#     Unit[] columns that are next to them.

.addUnitCols <- function( colIndices, SDRF ) {

    

    # Get the indices of unit columns. 

    unitCols <- unlist(

        

        # Go through each column index provided.

        sapply(

            colIndices,

            function( colNumber ) {

                    

                # Look at the column to the right of it.

                nextCol = colNumber + 1



                # Only try if this is not the very last column.

                if( nextCol <= ncol( SDRF ) ) {



                    # If it's a unit column, return the index.

                    if( grepl( "Unit", SDRF[ 1, nextCol ] ) ) {

                        nextCol

                    }

                }

            }

        )

    )



    # Combine the vector of unit column indices with the original one, and then

    # sort them so the unit column indices are next to the correct non-unit

    # column indices. E.g. if "Characteristics[ age ]" was the 2nd column out

    # of 5, that means its "Unit[ time unit ]" column is the 3rd column.

    # Combining the unit and non-unit column indices gives:

    #     1, 2, 4, 5, 3

    # So then we sort it so it becomes:

    #     1, 2, 3, 4, 5

    # Now when we use this vector of indices to select the Characteristics

    # columns we'll get them in the right order.

    allColIndices <- sort( c( colIndices, unitCols ) )



    return( allColIndices )

}





.mergeUnits <- function( subsetSDRF ) {



    # Find the unit columns.

    unitCols <- grep( "Unit", subsetSDRF[ 1, ] )

    

    if( length( unitCols ) > 0 ) {

        

        cat( "Merging unit columns...\n" )

    

        # Create some new merged columns.

        mergedCols <- data.frame( 

            sapply( 

                unitCols,

                function( unitCol ) {

                    paste( subsetSDRF[ , unitCol - 1 ], subsetSDRF[ , unitCol ] )

                }

            ),

            stringsAsFactors = FALSE

        )



        # Get the indices of the columns for which these unit columns apply.

        valueCols <- unitCols - 1

        

        # Replace the first row (header) with the one from the original SDRF.

        mergedCols[ 1 , ] <- subsetSDRF[ 1, valueCols ]

        

        # Replace the original value columns with the new merged columns.

        subsetSDRF[ , valueCols ] <- mergedCols



        # Delete the Unit columns.

        subsetSDRF <- subsetSDRF[ , -unitCols ]



        cat( "Finished merging unit columns.\n" )

    

        return( subsetSDRF )



    } else {



        # If there aren't any Unit columns, just return without doing anything.

        cat( "No Unit columns found.\n" )

        return( subsetSDRF )

    }

}

print(1)
experimentXMLlist <- parseAtlasConfig( opt$configuration )
experimentXMLlist
print(2)
atlasSDRF <- parseSDRF( opt$sdrf, experimentXMLlist$experimentType )
atlasSDRF
print(3)
allAnalytics <- experimentXMLlist$allAnalytics
allAnalytics
print(4)
expressionsFile <- file.path( opt$countstsv )
expressionsFile
print(5)
expressions <- read.delim( expressionsFile, header=TRUE, stringsAsFactors=FALSE, row.names = 1 )
summary(expressions)
print(6)

experimentSummary<-SimpleList( lapply( allAnalytics, function(analytics) {
  analyticsSDRF <- createAnalyticsSDRF( analytics, atlasSDRF )
  sumExp<-createSummarizedExperiment( expressions, analyticsSDRF, opt$methods )
  return(sumExp)
}) )
print("experiment summary 6.1")
head(experimentSummary)
study_batch<-as.factor(experimentSummary$rnaseq[[opt$batch]])
print(study_batch)
print("expe sum rnaseq 6.2")
experimentSummary$rnaseq
print("opt$batch 6.3")
opt$batch
print("study  batchv 7")
study_batch<-as.factor(experimentSummary$rnaseq[[opt$batch]])
study_batch
print("experimentSummary$rnaseq@metadata 8")
experimentSummary$rnaseq@metadata<-list(batch=study_batch)
experimentSummary$rnaseq@metadata
print("experimentSummary$rnaseq$batch")
experimentSummary$rnaseq$batch
print("experimentSummary$rnaseq$batch 9")

experimentSummary$rnaseq$batch<-study_batch
print("9.1")
experimentSummary$rnaseq$batch
print(10)
if( opt$remove_single_batch_genes ) {
  experimentSummary$rnaseq<-remove_genes_with_expression_on_single_batch(experimentSummary$rnaseq)
}


# Save the object to a file.
save( experimentSummary, file = opt$output )
