```{r}
feature_data |>
  mutate(feature_type = str_to_title(feature_type)) |>
  count(feature_type) |>
  arrange(desc(n)) |>  
  ggplot(aes(x = reorder(feature_type, -n),
             y = n)) +
  geom_bar(stat = "identity", 
           fill = "#1f618d",
           alpha = 0.8) +
  labs(title = "Frequency of Annotation Types",
       x = "Annotation Type",
       y = "Frequency") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```


```{r}
ggplot(length_data, 
       aes(x = log(length))) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 0.1, 
                 fill = "#1f618d",
                 color = "#1f618d",
                 alpha = 0.8) +
  labs(x = "Amino Acid Length",
       y = "Density") +
  theme_classic()
```

One might be concerned that by adjusting the length in this way, we would end up removing almost everything when an annotation spans a large area in the beginning of the sequence but extents beyond position 1024. Thus, we look into the adjusted length of the cropped sequences.

```{r}
feature_data |>
  mutate(width = end_position - start_position) |>
  filter(length > 1024) |>
  count(accession, length_adj) |>
  ggplot(aes(x = length_adj)) +
  geom_boxplot(fill = "#1f618d",
               alpha = 0.8) + 
  theme_classic()
```

It seems fair.

```{r}
feature_data |>
  mutate(width = end_position - start_position) |>
  filter(length > 1024) |>
  count(accession, length_adj) |>
  select(length_adj) |>
  pull() |>
  quantile(probs = c(0.1, 0.25, 0.5, 0.75, 0.90, 0.95))
```

For the sequences longer than 1024 corresponding to 2286 sequences, 25% will have an adjusted length of 854 or longer. Before filtering away the annotations, which will not be considered, we will take a look at what is removed.

```{r}
feature_data |>
  filter(!is.na(start_position)) |>
  mutate(keep = case_when(start_position <= length_adj ~ "Yes",
                          start_position > length_adj ~ "No")) |>
  mutate(feature_type = str_to_title(feature_type)) |>
  ggplot(aes(x = feature_type,
             fill = keep)) +
  geom_bar(position = "fill",
           alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Feature type",
       y = "Proportion",
       fill = "Keep") + 
  scale_fill_manual(values = c("#990000", "#1f618d"))
```


## Overview of annotations for smaller data set (100 sequence)

We will look at a subset of the data used for initial studies. First, let's load the accesion numbers.

```{r}
accessions_sp_100 <- read.table(file = "../../data/100_small/additional/UP000005640_9606_sp_100_acc.txt") |>
  rename(accession = V1)
```

We will filter out the sequences.

```{r}
feature_data_sp_100_filtered <- feature_data_sp_filtered |>
  filter(accession %in% accessions_sp_100$accession)
```

The annotations for the 100 sequences are saved.

```{r}
saveRDS(feature_data_sp_100_filtered, file = "../../data/100_small/annotations/features_sp_100_filtered.rds")
```

We take a look at which annotations are not included:

```{r}
unique(feature_data_sp_filtered$feature_type) %in% annotation_types$feature_type
```

All annotations after filtering are represented among the 100 sequences. We take a look at the length of the sequences.

```{r}
length_data_sp_100 <- length_data_sp |>
  filter(accession %in% accessions_sp_100$accession)
```

We count the number of sequence larger than 1024.

```{r}
n <- sum(length_data_sp_100$length > 1024)

print(paste("Sequences larger than 1024:", n))

"#A2B5CD" "#003366" "#1f618d"
```















