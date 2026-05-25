# Load the jsonlite package
#install.packages("jsonlite")
#install.packages("purr")
library(jsonlite)
library(purrr)
library(dplyr)
# Read the JSON file
json_data <- fromJSON("sample.stripy.json")

file_name <- basename("sample.stripy.json")
wgs_id <- sub("\\.[^.]*$", "", sub("\\.[^.]*$", "",sub("\\.[^.]*$", "",sub("\\.[^.]*$", "", file_name))))

# Remove the elements "Checksum" and "JobDetails" from the list
json_data <- json_data[!names(json_data) %in% c("Checksum", "JobDetails")]

# Extract the sub-elements
sub_elements <- unlist(json_data, recursive = FALSE)

# Remove the prefix "GenotypingResult." from all sub-elements
sub_elements_names <- gsub("^GenotypingResults\\.", "", names(sub_elements))

# Rename the sub-elements
names(sub_elements) <- sub_elements_names

# Convert the sub-elements back to JSON format
modified_json <- toJSON(sub_elements)

# Write the modified JSON data to a file
writeLines(modified_json, "modified_data.json")

json_data <- fromJSON("modified_data.json")

# Iterate through each element and remove the "SVG" attribute
for (element_name in names(json_data)) {
  json_data[[element_name]]$SVG <- NULL
}

# Convert the modified JSON data back to JSON format
modified_json <- toJSON(json_data)

# Write the modified JSON data to a file
writeLines(modified_json, "modified_data_no_svg.json")

json_data <- fromJSON("modified_data_no_svg.json")

sub_elements <- unlist(json_data, recursive = FALSE)
sub_elements <- unlist(sub_elements, recursive = FALSE)

sub_elements <- Filter(function(x) length(x) > 0, sub_elements)

# Remove NA values from all attributes in sub_elements
for (attr in names(sub_elements)) {
  if (is.numeric(sub_elements[[attr]])) {
    # For numeric vectors, remove NA values
    sub_elements[[attr]] <- sub_elements[[attr]][!is.na(sub_elements[[attr]])]
  } else if (is.character(sub_elements[[attr]])) {
    # For character vectors, remove NA values and empty strings
    sub_elements[[attr]] <- sub_elements[[attr]][sub_elements[[attr]] != "" & !is.na(sub_elements[[attr]])]
  }
}

# Filter out elements with NA values
sub_elements <- sub_elements[sapply(sub_elements, function(x) !any(is.na(unlist(x))))]

final_1 <- sub_elements

#sub_elements <- final_1

sub_elements <- unlist(sub_elements, recursive = FALSE)

sub_elements <- Filter(function(x) length(x) > 0, sub_elements)

sub_elements <- sub_elements[sapply(sub_elements, function(x) !any(is.na(unlist(x))))]

sub_elements <- unlist(sub_elements, recursive = FALSE)

# Function to split elements with multiple values such as CI
split_elements <- function(x) {
  if (length(x) > 1) {
    suffix <- 1:length(x)
    names(x) <- paste0(names(x), suffix)
    return(as.list(x))
  } else {
    return(as.list(x))
  }
}

# Apply the function to split elements
sub_elements <- purrr::map(sub_elements, split_elements)

# Flatten the list
sub_elements <- unlist(sub_elements, recursive = FALSE)

# Remove the Allele order from all sub-elements and Metadata.TotalOf and the digits at the end
sub_elements_names <- gsub("\\.Alleles\\d+\\.", ".Alleles.", names(sub_elements))
names(sub_elements) <- sub_elements_names

sub_elements_names <- gsub("Metadata\\.TotalOf(\\w*)\\d*", "\\1", names(sub_elements))
names(sub_elements) <- sub_elements_names

sub_elements_names <- gsub("\\Reads\\d*", ".Reads", names(sub_elements))
names(sub_elements) <- sub_elements_names

sub_elements_names <- gsub("(?<=Reads)\\d*$", "", names(sub_elements), perl = TRUE)
names(sub_elements) <- sub_elements_names

sub_elements_names <- gsub("(?<=Filter)\\d*$", "", names(sub_elements), perl = TRUE)
names(sub_elements) <- sub_elements_names

sub_elements_names <- gsub("(?<=Flag)\\d*$", "", names(sub_elements), perl = TRUE)
names(sub_elements) <- sub_elements_names

# Create a new JSON object with sample_name as the first key-value pair
new_sub_elements <- list(wgs_id = wgs_id)
new_sub_elements <- c(new_sub_elements, sub_elements)

###########################Numbering Coloumns for Sex and Non-sex Alleles##########

temp_new_sub_elements <- new_sub_elements

new_sub_elements <- temp_new_sub_elements

# Identify duplicate column names
table_names <- table(names(new_sub_elements))
duplicate_names <- names(table_names)[table_names > 1]

# Create a counter to keep track of the appended number for each unique column name
counter <- rep(1, length(duplicate_names))

# Loop through the list and rename duplicate columns
for (i in seq_along(new_sub_elements)) {
  # Get the current column name
  current_name <- names(new_sub_elements)[i]
  if (any(duplicate_names == current_name)) {
    # Get the index of the duplicate column name in the list of duplicates
    idx <- which(duplicate_names == current_name)
    # Append the counter to the current column name
    new_name <- paste0(current_name, counter[idx])
    # Increment the counter for the current column name
    counter[idx] <- counter[idx] + 1
    # Update the column name in the list
    names(new_sub_elements)[i] <- new_name
  }
}

#####################################################################

# Convert JSON data to a data frame
df <- as.data.frame(new_sub_elements)

# Identify the columns containing "Max" and "Min"
max_cols <- grep("\\.Max(\\.|\\d*$)", names(df), value = TRUE)
min_cols <- grep("\\.Min(\\.|\\d*$)", names(df), value = TRUE)

# Create new columns by merging Max and Min columns
for (i in seq_along(min_cols)) {
  min_col <- min_cols[i]
  max_col <- sub("\\.Min", ".Max", min_col)
  
  # Extract common column name
  common_name <- sub("\\.Min(\\.\\d+)?$", "", min_col)
  
  # Extract numeric suffix from max_col
  suffix <- sub("^.*\\.Max(\\.\\d+)?$", "\\1", max_col)
  suffix <- gsub("\\.", "", suffix)  # Remove the dot
  # Check if the suffix is empty (indicating only one pair of Max and Min)
  if (suffix == "") {
    # Create new column with common name and suffix
    new_col <- paste0(common_name, ":", 1)
    # Merge the values of Min and Max columns and assign to the new column
    df[[new_col]] <- paste(df[[min_col]], df[[max_col]], sep = "-")
  } else {
    # Extract the suffixes
    suffixes <- unique(gsub("\\.Max", "", suffix))
    # Loop over the suffixes and create new columns for each pair of Min and Max
    for (s in suffixes) {
      # Create new column with common name and suffix
      new_col_suffix <- paste0(common_name, ":", s)
      # Merge the values of Min and Max columns with the corresponding suffix and assign to the new column
      df[[new_col_suffix]] <- paste(df[[min_col]], df[[max_col]], sep = "-")
    }
  }
}

# Remove the original Max and Min columns
df <- df[, !names(df) %in% c(min_cols, max_cols), drop = FALSE]

# Sort the columns alphabetically, keeping the first column as "Samplename"
df <- df[, c("wgs_id", sort(names(df)[-1]))]

# Add "graph_hg38_" prefix to all column names except the first one
names(df)[-1] <- paste0("graph_hg38.", names(df)[-1])

#Add number on ChrX loci in Males

# Define the patterns to match
patterns <- c(
  "Alleles.IsPopulationOutlier",
  "Alleles.PopulationZscore",
  "Alleles.Range",
  "Alleles.Repeats",
  "Flanking.Reads",
  "Inrepeat.Reads",
  "Spanning.Reads"
)

# Get the names of the columns in the DataFrame
col_names <- names(df)

# Iterate through the patterns
for (pattern in patterns) {
  # Find the columns that match the pattern
  matching_cols <- grep(pattern, col_names, value = TRUE)
  
  # Append "1" to the matching column names if there's no digit already present
  for (col in matching_cols) {
    if (!grepl("\\d$", col)) {
      new_col <- paste0(col, "1")
      names(df)[names(df) == col] <- new_col
    }
  }
}

# Write the flattened data frame to a CSV file
write.csv(df, file = "sample.stripy.csv", row.names = FALSE)
















