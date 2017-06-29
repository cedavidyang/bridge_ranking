import numpy as np

def compute_cost(bridge_db, total_delay, total_distance, t):
    """compute risk according to Saydam and Frangopol (2011)"""
    TRUCK2TRAFFIC = 0.12    # Average Daily Truck Traffic (ADTT) to Average Daily Traffic (ADT)
    TRUCK_COMP = 26.97    # average compensation for truck drivers, USD/h
    CAR_OCCUPY = 1.5        # average vehicle occupancies for cars
    TRUCK_OCCUPY = 1.05     # average vehicle occupancies for trucks
    CAR_WAGE = 22.82        # average wage of car drivers, USD/h
    DISCOUNT_RATE = 0.00    # discount rate
    REBUILD_COST = 894      # rebuild cost, USD/m2
    CAR_RUN_COST = 0.08     # running cost for cars, USD/km
    TRUCK_RUN_COST = 0.375  # running cost for trucks, USD/km
    CARGO_VALUE = 4.        # time value of a cargo

    total_cost = 0.
    for (name, lat, long, length, width, deck_cs0, super_cs0, sub_cs0, detour, onlink) in bridge_db:
        cost_reb = REBUILD_COST * length*width * (1+DISCOUNT_RATE)**t
        cost_time = (total_delay/3600.*(1.-TRUCK2TRAFFIC)*CAR_OCCUPY*CAR_WAGE +
            total_delay/3600.*TRUCK2TRAFFIC*(TRUCK_COMP*TRUCK_OCCUPY+CARGO_VALUE))*(1+DISCOUNT_RATE)**t
        cost_run = (total_distance/1000.*(1.-TRUCK2TRAFFIC)*CAR_RUN_COST +
                total_distance/1000.*TRUCK2TRAFFIC*TRUCK_RUN_COST)*(1+DISCOUNT_RATE)**t
        total_cost += cost_reb+cost_time+cost_run

    return total_cost

def bridge_cost(bridge_db, t):
    """compute bridge cost according to Saydam and Frangopol (2011)"""
    DISCOUNT_RATE = 0.00    # discount rate
    REBUILD_COST = 894      # rebuild cost, USD/m2
    total_cost = 0.
    for (name, lat, long, length, width, deck_cs0, super_cs0, sub_cs0, detour, onlink) in bridge_db:
        cost_reb = REBUILD_COST * length*width * (1+DISCOUNT_RATE)**t
        total_cost += cost_reb

    return total_cost


def social_cost(total_delay, total_distance, t):
    """compute social cost according to Saydam and Frangopol (2011)"""
    TRUCK2TRAFFIC = 0.12    # Average Daily Truck Traffic (ADTT) to Average Daily Traffic (ADT)
    TRUCK_COMP = 26.97    # average compensation for truck drivers, USD/h
    CAR_OCCUPY = 1.5        # average vehicle occupancies for cars
    TRUCK_OCCUPY = 1.05     # average vehicle occupancies for trucks
    CAR_WAGE = 22.82        # average wage of car drivers, USD/h
    DISCOUNT_RATE = 0.00    # discount rate
    CAR_RUN_COST = 0.08     # running cost for cars, USD/km
    TRUCK_RUN_COST = 0.375  # running cost for trucks, USD/km
    CARGO_VALUE = 4.        # time value of a cargo

    cost_time = (total_delay/3600.*(1.-TRUCK2TRAFFIC)*CAR_OCCUPY*CAR_WAGE +
        total_delay/3600.*TRUCK2TRAFFIC*(TRUCK_COMP*TRUCK_OCCUPY+CARGO_VALUE))*(1+DISCOUNT_RATE)**t
    cost_run = (total_distance/1000.*(1.-TRUCK2TRAFFIC)*CAR_RUN_COST +
            total_distance/1000.*TRUCK2TRAFFIC*TRUCK_RUN_COST)*(1+DISCOUNT_RATE)**t
    total_cost = cost_time+cost_run

    return total_cost
