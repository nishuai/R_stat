###catch the best bet
odds=c(1.61, 3.9, 6.5)
sum(1/odds)

bets=c(500, 40, 58)

0.61* (4+6)
final_bet=function(bets, odds){
    rates=1/odds
    safe_deduction=rates*(min(bets/rates))
    bets-safe_deduction
}

bet2=final_bet(bets,odds)


sum(bets)-bets*odds
sum(bet2)-bet2*odds
sum(bet2)-sum(bets)
