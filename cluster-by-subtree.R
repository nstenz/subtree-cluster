#!/usr/bin/Rscript
library(ape);
library(phangorn);

prep_tree = function(target_tree) {

	# Read tree in
	tree = read.tree(target_tree);

	tree = midpoint(tree);
	tree = unroot(tree);
		
	# Plot tree
	plot(tree, type="phylogram", cex=0.2, main="Initial Tree, Midpoint Rooted, then Unrooted");

	return(tree);
}

reduce_subtrees = function(subtrees, species) {

	# Account for missing subtrees due to subtrees() not treating tree as unrooted
	missed_trees = vector("list", );
	class(missed_trees) = "multiPhylo";
	for (subtree in subtrees) {

		# Skip full tree
		#if (all.equal.phylo(subtree, tree)) {
		if (all.equal.phylo(subtree, working_tree)) {
			next;
		}

		# Remove this subtree's tips from the full tree
		#tip_tree = drop.tip(tree, subtree$tip.label);
		if (length(intersect(working_tree$tip.label, subtree$tip.label)) != length(working_tree$tip.label) - 1) {
			tip_tree = drop.tip(working_tree, subtree$tip.label);
		}
		else {
			next;
		}

		# Check if this subtree exists in the current subtree list
		in_current_subtrees = 0;
		for (subtree in subtrees) {
			if (all.equal.phylo(tip_tree, subtree)) {
				in_current_subtrees = in_current_subtrees + 1;  
				break;
			}
		}   

		# Add to list of missed trees if we don't have this subtree already
		if (!in_current_subtrees) {
			missed_trees[[length(missed_trees) + 1]] = tip_tree;        
		}   
	}
	#subtrees = c(subtrees, missed_trees);

	# Remove trees which don't include all species
	#cat(length(subtrees), "initial subtrees\n");
	for (index in length(subtrees):1) {
		subtree = subtrees[[index]];
		subtree_species = gsub("\\|.*$", "", subtree$tip.label, perl=T);
		subtree_species = unique(subtree_species);

		# Subtree does not include all species, remove
		if (length(subtree_species) != length(species)) {
			subtrees[[index]] = NULL;
		}
	}

	# Calculate smallest subtrees which don't include another subtree
	subtrees_reduced = subtrees;
	for (index in length(subtrees_reduced):1) {
		
		for (index2 in 1:length(subtrees_reduced)) {
			# Skip the same tree
			if (index2 == index) {
				next;
			}
		
			# Extract tip labels of subtrees
			subtree1_tip_labels = subtrees_reduced[[index]]$tip.label;
			subtree2_tip_labels = subtrees_reduced[[index2]]$tip.label;
			
			# Find shared tips
			shared_tips = intersect(subtree1_tip_labels, subtree2_tip_labels);

			# Remove the larger tree if there is overlap
			if (length(shared_tips) == length(subtree1_tip_labels)) {
				subtrees_reduced[[index2]] = NULL;
				break;
			}
			else if (length(shared_tips) == length(subtree2_tip_labels)) {
				subtrees_reduced[[index]] = NULL;
				break;
			}

#			# Remove the larger tree if there is overlap
#			if (length(shared_tips) > 0 && length(subtree1_tip_labels) > length(subtree2_tip_labels)) {
#				subtrees_reduced[[index]] = NULL;
#				break;
#			} else if (length(shared_tips) > 0 && length(subtree1_tip_labels) < length(subtree2_tip_labels)) {
#				subtrees_reduced[[index2]] = NULL;
#				break;
#			}
#			else if (length(shared_tips) > 0) {
#				cat("i think this should never happen\n");
#				q();
#			}
		}
	}

	# Enlarge subtrees
	for (index in 1:length(subtrees)) {
		
		contained_subtree_count = 0;
		contained_subtree_index = 1L;

		#for (subtree_reduced in subtrees_reduced) {
		for (index2 in 1:length(subtrees_reduced)) {

			# Extract tip labels of subtrees
			subtree1_tip_labels = subtrees[[index]]$tip.label;
			subtree2_tip_labels = subtrees_reduced[[index2]]$tip.label;
			
			# Find shared tips
			shared_tips = intersect(subtree1_tip_labels, subtree2_tip_labels);

			# Same tree, skip
			if (length(shared_tips) == length(subtree1_tip_labels) && length(shared_tips) == length(subtree1_tip_labels)) {
				next;
			}

			# Tips are shared between the two trees
			if (length(shared_tips) > 0) {
				contained_subtree_index = index2;
				contained_subtree_count = contained_subtree_count + 1;
			}
		}

		# Expand the reduced subtree to the larger one if it only contained a single subtree
		if (contained_subtree_count == 1) {

			# Check which tips we've added
			added_tips = setdiff(subtrees[[index]]$tip.label, subtrees_reduced[[contained_subtree_index]]$tip.label);

			# Additional check for enlarging subtree
			if (!(length(grep("^(Grai)|(hirD)", added_tips, perl=T)) > 0 && length(grep("^(Garb)|(hirA)", added_tips, perl=T)))) {
				subtrees_reduced[contained_subtree_index] = subtrees[index];
			}
		}
	}

#	for (index in 1:length(subtrees_reduced)) {
#		if (length(subtrees_reduced[[index]]$tip.label) > 3) {
#			#subtrees_reduced[[index]] = unroot(midpoint(subtrees_reduced[[index]]));
#			subtrees_reduced[[index]] = midpoint(subtrees_reduced[[index]]);
#		}
#		else {
#			subtrees_reduced[[index]] = midpoint(subtrees_reduced[[index]]);
#		}
#	}

	return(subtrees_reduced);
}

# Parse command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
target_tree = args[1];
#target_gene = gsub("\\.ph", "", target_tree);
target_gene = gsub("(.*)\\..*$", "\\1", target_tree);

# Check that user gave a tree file
if (!exists("target_tree")) {
	cat("You must specify a tree file as input.\n");
	q("no", status = 1);
} else if (!file.exists(target_tree)) {
	cat(paste0("Could not locate '", target_tree, "' perhaps you made a typo?\n"));
	q("no", status = 1);
}

# Open output pdf
#pdf(paste0(target_gene, ".pdf"));
pdf(paste0(target_gene, ".pdf"), width=10, height=10);
tree = prep_tree(target_tree);

#tip_colors = rep("black", length(tree$tip.label));
tip_colors = rep(1, length(tree$tip.label));

# Vector of trees that need to be parsed
unparsed_trees = vector("list", );
class(unparsed_trees) = "multiPhylo";
unparsed_trees[[1]] = tree;

#start = 1;
clade_color = 1;
subtree_index = 1;
subtree_count = 1;
while (subtree_index <= length(unparsed_trees)) {

	# Tree we are currently working on
	working_tree = unparsed_trees[[subtree_index]];
	plot(working_tree, type="phylogram", main=c("working tree"), cex=0.2);

	# Extract names of all species
	species = gsub("\\|.*$", "", working_tree$tip.label, perl=T);
	species = unique(species);

	# Extract all subtrees of our tree
	#subtrees = subtrees(working_tree, wait=T);
	subtrees = subtrees(working_tree);
	subtrees_reduced = reduce_subtrees(subtrees, species);

	# Output final subtrees
	tips_in_subtree = vector();
	for (index in 1:length(subtrees_reduced)) {
		subtree = subtrees_reduced[[index]];

		# Color code clades
		overlap = match(subtree$tip.label, tree$tip.label);
		tips_in_subtree = union(tree$tip.label[overlap], tips_in_subtree);
		for (match in overlap) {
			tip_colors[match] = clade_color;
		}

		# Remove node IDs which are showing up for some reason
		out_tree = subtree;
		out_tree$node.label = NULL;

		filename = paste0(target_gene, "_clusters.tre");
		if (clade_color == 1) {
			write.tree(out_tree, file=filename);
		}
		else {
			write.tree(out_tree, file=filename, append=T);
		}

		# Increment clade color
		clade_color = clade_color + 1;
		#plot(subtree, type="phylogram", main=c("subtree", index));
	}

	# Determine which tips are not in a subtree
	tips_not_in_subtree = setdiff(working_tree$tip.label, tips_in_subtree);
	tips_not_in_subtree_indices = match(tips_not_in_subtree, tree$tip.label);

	# Split tips which aren't in a subtree with all members into their own respective subtrees

	tree_index = 1;
	tips_not_in_subtree_trees = vector("list", );
	class(tips_not_in_subtree_trees) = "multiPhylo";
	start_tip = tips_not_in_subtree_indices[1];

	if (length(tips_not_in_subtree_indices > 1)) {

		# Group tips by their ID
		tip_ranges = vector("character", );
		for (i in 1:(length(tips_not_in_subtree_indices) - 1)) {
			index = tips_not_in_subtree_indices[i];
			next_index = tips_not_in_subtree_indices[i + 1];

			#cat("index", index, "next_index", next_index, "\n");

			# Clade ends if indices aren't 1 apart or if we are at the last index in the list
			if (next_index != index + 1 || i == length(tips_not_in_subtree_indices) - 1) {

				# Prevents skipping of last tip
				if (i == length(tips_not_in_subtree_indices) - 1) {
					index = next_index;
				}
			
				# Set range of tips
				tip_ranges = c(tip_ranges, paste0(start_tip, ":", index));	
				start_tip = next_index;
			}
		}

		# Check if the tips in the range are monophyletic
		for (tip_range in tip_ranges) {
			range = unlist(strsplit(tip_range, ":"));	

			start = as.integer(range[[1]]);
			end = as.integer(range[[2]]);

			# Add to list of unparsed trees if monophyletic
			if (is.monophyletic(tree, tree$tip.label[start:end])) {

				# Get tips not in the subtree
				inverse_clade_tips = tree$tip.label[seq(-1 * start, -1 * end, -1)];

				# Extract clade by removing tips not in subtree
				tips_not_in_subtree_trees[[tree_index]] = drop.tip(working_tree, inverse_clade_tips);

				# Set file name and output tree
				filename = paste0(target_gene, "_missing", tree_index, ".tre");

				# Add this subtree to out list of subtrees to check
				unparsed_trees[[subtree_count + 1]] = tips_not_in_subtree_trees[[tree_index]];
				subtree_count = subtree_count + 1;

				# Iterate counters for next clade
				tree_index = tree_index + 1;
			}
			else {
				print(tip_range);

				# Determine which tips form monophyletic clades
				mono = vector("character", );
				for (i in start:(end - 1)) {
					tip1 = i;

					for (j in (i + 1):end) {
						tip2 = j;

						if (is.monophyletic(tree, c(tree$tip.label[tip1], tree$tip.label[tip2]))) {
							mono = c(mono, paste0(tip1, ":", tip2));	
							cat(tree$tip.label[tip1], tree$tip.label[tip2],"\n");
						}
					}
				}
				print(mono);
				#q();

				# Remove overlapping monophyletic clades if a larger one exists
				mono_reduced = mono;
				for (index in (length(mono_reduced)):1) {
					for (index2 in 1:(length(mono_reduced))) {
						if (index == index) {
							next;
						}

						range1 = unlist(strsplit(mono_reduced[index], ":"));
						range2 = unlist(strsplit(mono_reduced[index2], ":"));

						# TODO: REDUCE
						start1 = as.integer(range1[[1]]);
						end1 = as.integer(range1[[2]]);

						start2 = as.integer(range2[[1]]);
						end2 = as.integer(range2[[2]]);

						range1 = start1:end1;
						range2 = start2:end2;

						overlap = intersect(range1, range2);

						# Check if there is overlap
						if (overlap) {
							cat("overlapping monophyletic:", overlap,"\n");

							# Remove smaller 
							if (length(range2) < length(range1)) {
								mono_reduced[[index2]] = NULL;
							}
							else {
								mono_reduced[[index]] = NULL;
							}
							break;
						}
					}
				}

				# Handle single tips not in a monophyletic group
				for (index in start:end) {
					in_group = 0;
					for (group in mono_reduced) {
							range1 = unlist(strsplit(group, ":"));

							# TODO: REDUCE
							start1 = as.integer(range1[[1]]);
							end1 = as.integer(range1[[2]]);

							range1 = start1:end1;

						if (index %in% range1) {
							in_group = 1;				
						}
					}

					if (!in_group) {
						tip_colors[index] = clade_color;

						# Increment clade color
						clade_color = clade_color + 1;
					}
				}

				# Add monophyletic groups to unparsed subtrees
				for (group in mono_reduced) {
					range1 = unlist(strsplit(group, ":"));

					# TODO: REDUCE
					start1 = as.integer(range1[[1]]);
					end1 = as.integer(range1[[2]]);

					range1 = start1:end1;

					# Get tips not in the subtree
					inverse_clade_tips = tree$tip.label[seq(-1 * start1, -1 * end1, -1)];

					# Extract clade
					clade = drop.tip(working_tree, inverse_clade_tips);
					tips_not_in_subtree_trees[[tree_index]] = clade;

					# Set file name and output tree
					filename = paste0(target_gene, "_missing", tree_index, ".tre");

					# Add this subtree to out list of subtrees to check
					unparsed_trees[[subtree_count + 1]] = tips_not_in_subtree_trees[[tree_index]];
					subtree_count = subtree_count + 1;

					# Iterate counters for next clade
					tree_index = tree_index + 1;
				}
			}
		}
	}
	else {
		for (index in tips_not_in_subtree_indices) {
			tip_colors[index] = clade_color;
		}

		# Increment clade color
		clade_color = clade_color + 1;
	}

	subtree_index = subtree_index + 1;
}

# Set up clade-color coding
colors = grDevices::colors()[grep('(gr(a|e)y)|(black)', grDevices::colors(), invert = T)]
#colors = colorRampPalette(c("red", "orange", "yellow", "green", "blue", "violet"))(max(tip_colors))
colors = sample(colors, length(colors));

tip_colors = colors[tip_colors];

#plot(tree, cex = 0.2, type="radial", main=c("full with clades"), show.tip.label=F);
#tiplabels(tree$tip.label, col=tip_colors, cex=0.2, bg="black");

# Using a black background seems to be the only way to reliably get a readable plot
par(bg="black");
#plot(tree, cex = 0.5, type="fan", main=c("full with clades"), tip.color=tip_colors, lab4ut="axial", edge.color="white", adj=0);
#plot(tree, cex = 0.5, type="fan", main=c("full with clades"), tip.color=tip_colors, lab4ut="axial", adj=0);
#plot(tree, cex = 0.5, type="fan", main=c("full with clades"), tip.color=tip_colors, lab4ut="axial", label.offset=0.01);
#plot(tree, cex = 0.5, type="fan", main=c("full with clades"), tip.color=tip_colors, lab4ut="axial", edge.color="white", label.offset=0.01);
plot(midpoint(tree), cex = 0.4, type="fan", tip.color=tip_colors, lab4ut="axial", edge.color="white", label.offset=0.01);
title(main=c("Clades"), col.main="white");
#tiplabels(tip_colors, cex = 0.3);
#tiplabels(cex = 0.3);

dev.off();
