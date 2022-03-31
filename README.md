Quantifying and correcting bias due to outcome dependent self-reported weights in longitudinal study of weight loss interventions
===================

Background
---

In the past few decades, the increasing rate of obesity has become a worldwide issue, leading to the increasing incidence of several diseases and economic burden. To motivate participants to lose weight and maintain healthy behavior, financial incentive interventions, such as lottery-based and direct payment incentives, are utilized and have been shown to be successful in practice. The objective of this paper is to propose a novel framework of methods to analyze longitudinal self-reported weight loss data and study the effectiveness of financial interventions.


Data
---

The Keep It Off study, a three-arm randomized controlled trial (RCT) with 189 participants in the analysis cohort, with daily self-reported weights, examined if the participants’ weight loss maintenance can be improved by financial incentives. 

![image](https://user-images.githubusercontent.com/38872447/161109621-8645cea4-9d22-4eed-890e-4c8e5ccb678f.png)



Methods
---

To overcome the challenges due to the informative reporting process in the real-world longitudinal data, we proposed a framework of methods to quantify the evidence of missing not at random due to the outcome-dependent self-reporting mechanism, and to conduct bias correction using estimating equation derived from pairwise composite likelihood.

![image](https://user-images.githubusercontent.com/38872447/161109563-61652411-37f6-4bb6-b73c-b4d1abaa21fd.png)



Implementation
---

In Stage II, the pairwise construction of likelihood comes with the price of higher computational cost, as the algorithm involves computation of likelihood constructed by all pairs of patients within a site. To alleviate this limitation, we implemented an algorithm with R calling C, which is about 50 times faster than using the R programming language alone. The R code and C code to implement the method can be found in this repo. 




