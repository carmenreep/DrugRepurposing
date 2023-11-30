import os
from flask import Flask, render_template, url_for, request, session,redirect
from drugapp import app
from drugapp.forms import Symptoms, MIMPhenoType
from drugapp.Monarch import *
from drugapp.DGIdb import *
from drugapp.drugsimilarity import *
from drugapp.combine_graphs import *
#from drugapp.neo4jlib import *
from drugapp.rdf2vec import *
import datetime
import pandas as pd
from pathlib import Path

# Make sure current working directory is always the drugapp folder:
current_directory = os.getcwd()
# if drugapp folder is not in the cwd (when app is started), path_out should first go into drugapp folder
if 'drugapp' not in current_directory:
    path_out = './drugapp/templates'
else:
    # if drugapp folder is in the cwd, change cwd until drugapp folder is reached
    current_directory = os.getcwd()
    while not (current_directory[-7:] == 'drugapp'):
        os.chdir("..")
        current_directory = os.getcwd()
        # if drugapp not in the current directory, break the while loop
        if len(current_directory) < 5:
            print('the drugapp folder is not in the cwd')
            break
    os.chdir("..") # one more up to go out of drugapp
    path_out = './drugapp/templates'

# make sure files in 'templates' are shown that are not the templates for creating the app
# and make sure latest file is at the top
path_out_sorted = sorted(Path(path_out).iterdir(), key=os.path.getmtime,reverse=True)
try:
    folders_output = [os.path.basename(d) for d in path_out_sorted if not (os.path.basename(d).startswith('about') or 
                                                                 os.path.basename(d).startswith('home') or 
                                                                 os.path.basename(d).startswith('layout') or 
                                                                 os.path.basename(d).startswith('symptoms'))] 
except StopIteration:
    print('StopIteration')

@app.route("/", methods = ['GET', 'POST'])
@app.route("/home", methods = ['GET', 'POST'])
def config():
    
    # Make sure current working directory is always the drugapp folder:
    current_directory = os.getcwd()
    # if drugapp folder is not in the cwd (when app is started), path should first go into drugapp folder
    if 'drugapp' not in current_directory:
        path = './drugapp/data'
    else:
        # if drugapp folder is in the cwd, change cwd until drugapp folder is reached
        current_directory = os.getcwd()
        while not (current_directory[-7:] == 'drugapp'):
            os.chdir("..")
            current_directory = os.getcwd()
            # if drugapp not in the current directory, break the while loop
            if len(current_directory) < 5:
                print('the drugapp folder is not in the cwd')
                break
        os.chdir("..") # one more up to go out of drugapp
        path = './drugapp/data'

    path_sorted = sorted(Path(path).iterdir(), key=os.path.getmtime, reverse=True) # make sure latest folder is at the top
    folderlist = [os.path.basename(d) for d in path_sorted] # make sure folders in 'data' are shown, which are the diseases
    
    form = MIMPhenoType()
    if form.validate_on_submit(): # if a mim number is inserted
        input_number = form.MIMPhenoType.data
        run_monarch(input_number)
        symptoms_name_lst, date, diseasename = symptom_list_today()
        session['symptoms_name_lst'] = symptoms_name_lst
        return redirect(url_for('symptoms', date=date, diseasename=diseasename))
    return render_template('home.html', folders = folderlist, form=form, folders_output=folders_output)

@app.route("/symptoms", methods = ['GET', 'POST'])
@app.route("/symptoms/<date>/<diseasename>", methods = ['GET', 'POST'])
def symptoms(date='None', diseasename='None'):
    form=Symptoms()
    form.symptoms.choices = session['symptoms_name_lst']
    if form.validate_on_submit():
        input_symptom = form.symptoms.data
        run_monarch_symptom(input_symptom, date)
        run_dgidb(date)
        run_drugsimilarity()
        run_combine_graphs(date)
        print(input_symptom)
        rdf2vec_general(input_symptom, date, diseasename)
        return render_template('{}; {}.html'.format(diseasename, input_symptom), folders_output=folders_output)
    return render_template('symptoms.html', title='Symptoms', form=form, diseasename=diseasename, folders_output=folders_output)

@app.route("/symptoms_from_folder/<diseasename>", methods = ['GET', 'POST'])
def symptoms_from_folder(diseasename='None'): # for when a disease folder is selected
    symptoms_name_lst, date, input_number= symptom_list_specified_folder(input_folder=diseasename)
    session['symptoms_name_lst'] = symptoms_name_lst
    return redirect(url_for('symptoms', date=date, diseasename=diseasename, input_number=input_number))

@app.route("/predictions/<filename>", methods = ['GET', 'POST'])
def predictions(filename='None'):
    return render_template(filename, folders_output=folders_output)
    
@app.route("/about")
def about():
    return render_template('about.html', title='About', folders_output=folders_output)