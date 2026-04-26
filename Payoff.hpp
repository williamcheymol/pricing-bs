#ifndef PAYOFF_HPP
#define PAYOFF_HPP

#include <algorithm> 
#include <cmath> 

/**
 * @class Payoff
 * @brief Classe de base abstraite pour représenter le payoff d'une option.
 */
class Payoff { 
public:
    /**
     * @brief Destructeur virtuel par défaut.
     */
    
    virtual ~Payoff() {}

    /**
     * @brief Calcule la valeur du payoff à l'échéance.
     * @param S Le prix du sous-jacent  à l'échéance.
     * @return La valeur de l'option (max(S-K, 0) ou max(K-S, 0)).
     */
    virtual double operator()(double S) const = 0;

    /**
     * @brief Calcule la condition au bord quand le prix du sous-jacent vaut 0.
     * @param t Le temps actuel.
     * @param T La maturité de l'option.
     * @param r Le taux d'intérêt sans risque.
     * @return La valeur de l'option quand S=0.
     */
    virtual double price_zero(double t, double T, double r) const = 0;

    /**
     * @brief Calcule la condition au bord quand le prix du sous-jacent vaut L.
     * @param t Le temps actuel.
     * @param T La maturité de l'option.
     * @param r Le taux d'intérêt sans risque.
     * @return La valeur de l'option quand S=L.
     */
    virtual double price_limit(double t, double T, double r) const = 0;
};

/**
 * @class Call
 * @brief Représente une option d'achat .
 * * Hérite de la classe Payoff. Le payoff est donné par max(S - K, 0).
 */
class Call : public Payoff {
private:
    double K_; 

public:
    /**
     * @brief Constructeur du Call.
     * @param K Le Strike.
     */
    Call(double K) : K_(K) {}

    /**
     * @brief Calcule le payoff du Call : max(S - K, 0).
     * @param S Prix du sous-jacent.
     * @return La valeur  du Call.
     */    
    double operator()(double S) const  {
        return std::max(0.0, S - K_);
    }

    /**
     * @brief Condition limite inférieure pour un Call.
     * * Si S = 0, le Call vaut 0 .
     * @return 0.0
     */
    double price_zero(double , double , double ) const  {
        return 0.0;
    }

    /**
     * @brief Condition limite supérieure pour un Call.
     * * Si S vaut L, le Call vaut S = K*exp(-r(t-T)).
     * @param t Temps actuel.
     * @param T Maturité.
     * @param r Taux sans risque.
     * @return   K * exp(-r(t-T))
    */
    double price_limit(double t, double T, double r) const {
        return K_ * std::exp(-r * (t - T));
    }
};

/**
 * @class Put
 * @brief Représente une option de vente .
 * * Hérite de la classe Payoff. Le payoff est donné par max(K - S, 0).
 */
class Put : public Payoff {
private:
    double K_; 

public:
    /**
     * @brief Constructeur du Put.
     * @param K Le Strike.
     */
    Put(double K) : K_(K) {}

    /**
     * @brief Calcule le payoff du Put : max(K - S, 0).
     * @param S Prix du sous-jacent.
     * @return La valeur du Put.
     */    
    double operator()(double S) const  {
        return std::max(0.0, K_ - S);
    }

    /**
     * @brief Condition limite inférieure pour un Put.
     * * Quand S = 0, le Put vaut K * exp(-r(T-t)).
     * @param t Temps actuel.
     * @param T Maturité.
     * @param r Taux sans risque.
     * @return K_ * std::exp(-r * (T - t)).
     */    
    double price_zero(double t, double T, double r) const {
        return K_ * std::exp(-r * (T - t));
    }

    /**
     * @brief Condition limite supérieure pour un Put.
     * * Quand S vaut L, l'option de vente vaut 0.
     * @return 0.0
     */    
    double price_limit(double , double , double  ) const {
        return 0.0;
    }
};

#endif